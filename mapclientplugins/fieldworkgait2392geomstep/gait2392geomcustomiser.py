"""
Gait2392 customisation
"""
import os
import numpy as np
import copy

from gias2.fieldwork.field import geometric_field
from gias2.mesh import vtktools
from gias2.musculoskeletal import mocap_landmark_preprocess
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfit
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfitscaling
from gias2.musculoskeletal import osim

import opensim

SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
TEMPLATE_OSIM_PATH = os.path.join(SELF_DIR, 'data/gait2392_simbody.osim')
OSIM_FILENAME = 'gait2392_simbody.osim'
OSIM_BODY_NAME_MAP = {'pelvis': 'pelvis',
                      'femur': 'femur_l',
                      'tibiafibula': 'tibia_l',
                      }
PELVIS_SUBMESHES = ('RH', 'LH', 'sac')
PELVIS_SUBMESH_ELEMS = {'RH': range(0, 73),
                        'LH': range(73,146),
                        'sac': range(146, 260),
                        }
PELVIS_BASISTYPES = {'tri10':'simplex_L3_L3','quad44':'quad_L3_L3'}
TIBFIB_SUBMESHES = ('tibia', 'fibula')
TIBFIB_SUBMESH_ELEMS = {'tibia': range(0, 46),
                        'fibula': range(46,88),
                        }
TIBFIB_BASISTYPES = {'tri10':'simplex_L3_L3','quad44':'quad_L3_L3'}

GEOM_DIR = 'geom/'
SACRUM_FILENAME = 'sacrum.vtp'
HEMIPELVIS_RIGHT_FILENAME = 'pelvis.vtp'
HEMIPELVIS_LEFT_FILENAME = 'l_pelvis.vtp'
FEMUR_LEFT_FILENAME = 'l_femur.vtp'
TIBFIB_LEFT_FILENAME = 'l_tibia.vtp'

class Gait2392GeomCustomiser(object):
    gfield_disc = (8,8)

    def __init__(self, config):
        self.config = config
        self._pelvisRigid = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self._hipRot = np.array([0.0, 0.0, 0.0])
        self._kneeRot = np.array([0.0, 0.0, 0.0])
        self.uniformScaling = 1.0
        self.pelvisScaling = 1.0
        self.femurScaling = 1.0
        self.petallaScaling = 1.0
        self.tibfibScaling = 1.0
        self.LL = None  # lowerlimb object
        self.osimmodel = None  # opensim model

        self._init_osim_model()

    def _init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)

    def _splitPelvisGFs(self):
        """
        Given a flattened pelvis model, create left hemi, sacrum,
        and right hemi meshes
        """
        gf = self.LL.models['pelvis'].gf
        lhgf = gf.makeGFFromElements(
                    'hemipelvis-left',
                    PELVIS_SUBMESH_ELEMS['LH'],
                    PELVIS_BASISTYPES
                    )
        sacgf = gf.makeGFFromElements(
                    'sacrum',
                    PELVIS_SUBMESH_ELEMS['sac'],
                    PELVIS_BASISTYPES
                    )
        rhgf = gf.makeGFFromElements(
                    'hemipelvis-right',
                    PELVIS_SUBMESH_ELEMS['RH'],
                    PELVIS_BASISTYPES
                    )
        return lhgf, sacgf, rhgf

    def _check_geom_path(self):
        """
        Check that the directory for geom meshes exists. If not, create it.
        """
        geom_dir = os.path.join(self.config['osim_output_dir'], GEOM_DIR)
        if not os.path.isdir(geom_dir):
            os.mkdir(geom_dir)
        
    def set_model_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """

        self.LL = bonemodels.LowerLimbLeftAtlas('left lower limb')
        self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
        self.LL.set_bone_gfield('femur', gfieldsdict['femur'])
        self.LL.set_bone_gfield('patella', gfieldsdict['patella'])
        self.LL.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula'])

    @property
    def pelvisRigid(self):
        return self._pelvisRigid

    @pelvisRigid.setter
    def pelvisRigid(self, value):
        if len(value)!=6:
            raise ValueError('input pelvisRigid vector not of length 6')
        else:
            self._pelvisRigid = np.array([value[0], value[1], value[2],
                                          _trimAngle(value[3]),
                                          _trimAngle(value[4]),
                                          _trimAngle(value[5]),
                                         ])

    @property
    def hipRot(self):
        return self._hipRot

    @hipRot.setter
    def hipRot(self, value):
        if len(value)!=3:
            raise ValueError('input hipRot vector not of length 3')
        else:
            self._hipRot = np.array([_trimAngle(v) for v in value])

    @property
    def kneeRot(self):
        if self.kneeDOF:
            return self._kneeRot[[0,2]]
        else:
            return self._kneeRot[[0]]

    @kneeRot.setter
    def kneeRot(self, value):
        if self.kneeDOF:
            self._kneeRot[0] = _trimAngle(value[0])
            self._kneeRot[2] = _trimAngle(value[1])
        else:
            self._kneeRot[0] = _trimAngle(value[0])

    @property
    def uniformScalingX(self):
        self._uniformScalingX = np.hstack([
                                self.uniformScaling,
                                self.pelvisRigid,
                                self.hipRot,
                                self.kneeRot
                                ])
        return self._uniformScalingX

    @uniformScalingX.setter
    def uniformScalingX(self, value):
        print value
        a = 1
        self._uniformScalingX = value
        self.uniformScaling = value[0]
        self.pelvisRigid = value[1]
        self.hipRot = value[2]
        self.kneeRot = value[3]

        # propagate isotropic scaling to each bone
        self.pelvisScaling = self.uniformScaling
        self.femurScaling = self.uniformScaling
        self.patellaScaling = self.uniformScaling
        self.tibfibScaling = self.uniformScaling

        self.lastTransformSet = self.uniformScalingX

    @property
    def perBoneScalingX(self):
        self._perBoneScalingX = np.hstack([
                                self.pelvisScaling,
                                self.femurScaling,
                                self.patellaScaling,
                                self.tibfibScaling,
                                self.pelvisRigid,
                                self.hipRot,
                                self.kneeRot
                                ])
        return self._perBoneScalingX

    @perBoneScalingX.setter
    def perBoneScalingX(self, value):
        a = 4
        self._perBoneScalingX = value
        self.pelvisScaling = value[0][1][0]
        self.femurScaling = value[0][1][1]
        self.patellaScaling = value[0][1][2]
        self.tibfibScaling = value[0][1][3]
        self.pelvisRigid = value[1]
        self.hipRot = value[2]
        self.kneeRot = value[3]
        self.lastTransformSet = self.perBoneScalingX

    def _save_vtp(self, gf, filename, bodycoordmapper):
        v, f = gf.triangulate(self.gfield_disc)
        v_local = bodycoordmapper(v)
        vtkwriter = vtktools.Writer(
                        v=v_local,
                        f=f,
                        filename=filename,
                        )
        vtkwriter.writeVTP()

    def cust_osim_pelvis(self):
        osim_pelvis = self.osimmodel.getBody(OSIM_BODY_NAME_MAP['pelvis'])
        pelvis = self.LL.models['pelvis']


        # update mass and inertial

        # update coordinate defaults 

        ## pelvis_tilt

        ## pelvis_list

        ## pelvis_rotation

        # update mesh
        lhgf, sacgf, rhgf = self._splitPelvisGFs()
        self._check_geom_path()

        ## sacrum.vtp
        sac_vtp_fn = os.path.join(self.config['osim_output_dir'], GEOM_DIR, SACRUM_FILENAME)
        self._save_vtp(sacgf, sac_vtp_fn, pelvis.acs.map_local)
        osim_pelvis.setDisplayGeometryFileName(sac_vtp_fn, SACRUM_FILENAME)

        ## pelvis.vtp
        rh_vtp_fn = os.path.join(self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_RIGHT_FILENAME)
        self._save_vtp(rhgf, rh_vtp_fn, pelvis.acs.map_local)
        osim_pelvis.setDisplayGeometryFileName(rh_vtp_fn, HEMIPELVIS_RIGHT_FILENAME)

        ## l_pelvis.vtp
        lh_vtp_fn = os.path.join(self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_LEFT_FILENAME)
        self._save_vtp(lhgf, lh_vtp_fn, pelvis.acs.map_local)
        osim_pelvis.setDisplayGeometryFileName(lh_vtp_fn, HEMIPELVIS_LEFT_FILENAME)


    def cust_osim_femur_left(self):
        osim_femur = self.osimmodel.getBody(OSIM_BODY_NAME_MAP['femur'])
        femur = self.LL.models['femur']

        # update hip_l joint 

        ## location in parent

        ## location in self

        # update coordinate defaults

        ## hip_flexion_l

        ## hip_adduction_l

        ## hip_rotation_l

        # update mesh l_femur.vtp
        self._check_geom_path()
        femur_vtp_fn = os.path.join(self.config['osim_output_dir'], GEOM_DIR, FEMUR_LEFT_FILENAME)
        self._save_vtp(femur.gf, femur_vtp_fn, femur.acs.map_local)
        osim_femur.setDisplayGeometryFileName(femur_vtp_fn, FEMUR_LEFT_FILENAME)

    def cust_osim_tibiafibula_left(self):
        osim_tibfib = self.osimmodel.getBody(OSIM_BODY_NAME_MAP['tibiafibula'])
        tibfib = self.LL.models['tibiafibula']
        # update knee_l joint

        ## location in parent

        ## location in self

        # update mesh l_tibia.vtp
        self._check_geom_path()
        tibfib_vtp_fn = os.path.join(self.config['osim_output_dir'], GEOM_DIR, TIBFIB_LEFT_FILENAME)
        self._save_vtp(tibfib.gf, tibfib_vtp_fn, tibfib.acs.map_local)
        osim_tibfib.setDisplayGeometryFileName(tibfib_vtp_fn, TIBFIB_LEFT_FILENAME)

    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(self.config['osim_output_dir'], TEMPLATE_OSIM_FILENAME)
            )