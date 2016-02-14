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

#=============================================================================#
SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
TEMPLATE_OSIM_PATH = os.path.join(SELF_DIR, 'data/gait2392_simbody.osim')
OSIM_FILENAME = 'gait2392_simbody.osim'
OSIM_BODY_NAME_MAP = {'pelvis': 'pelvis',
                      'femur-left': 'femur_l',
                      'femur-right': 'femur',
                      'tibiafibula-left': 'tibia_l',
                      'tibiafibula-right': 'tibia',
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
TIBIA_LEFT_FILENAME = 'l_tibia.vtp'
FIBULA_LEFT_FILENAME = 'l_fibula.vtp'

#=============================================================================#
# Opensim coordinate systems for bodies

# def update_femur_opensim_acs(femur_model):
#     femur_model.acs.update(
#         *bonemodels.model_alignment.createFemurACSOpenSim(
#             femur_model.landmarks['femur-HC'],
#             femur_model.landmarks['femur-MEC'],
#             femur_model.landmarks['femur-LEC'],
#             side=femur_model.side
#             )
#         )
# bonemodels.FemurModel.update_acs = update_femur_opensim_acs

# def update_tibiafibula_opensim_acs(tibiafibula_model):
#     tibiafibula_model.acs.update(
#         *bonemodels.model_alignment.createTibiaFibulaACSOpenSim(
#             tibiafibula_model.landmarks['tibiafibula-MM'],
#             tibiafibula_model.landmarks['tibiafibula-LM'],
#             tibiafibula_model.landmarks['tibiafibula-MC'],
#             tibiafibula_model.landmarks['tibiafibula-LC'],
#             side=tibiafibula_model.side
#             )
#         )
# bonemodels.TibiaFibulaModel.update_acs = update_tibiafibula_opensim_acs

def _splitTibiaFibulaGFs(tibfibGField):
    tib = tibfibGField.makeGFFromElements(
            'tibia',
            TIBFIB_SUBMESH_ELEMS['tibia'],
            TIBFIB_BASISTYPES,
            )
    fib = tibfibGField.makeGFFromElements(
            'fibula',
            TIBFIB_SUBMESH_ELEMS['fibula'],
            TIBFIB_BASISTYPES,
            )

    return tib, fib

def _splitPelvisGFs(pelvisGField):
    """
    Given a flattened pelvis model, create left hemi, sacrum,
    and right hemi meshes
    """
    lhgf = pelvisGField.makeGFFromElements(
                'hemipelvis-left',
                PELVIS_SUBMESH_ELEMS['LH'],
                PELVIS_BASISTYPES
                )
    sacgf = pelvisGField.makeGFFromElements(
                'sacrum',
                PELVIS_SUBMESH_ELEMS['sac'],
                PELVIS_BASISTYPES
                )
    rhgf = pelvisGField.makeGFFromElements(
                'hemipelvis-right',
                PELVIS_SUBMESH_ELEMS['RH'],
                PELVIS_BASISTYPES
                )
    return lhgf, sacgf, rhgf

#=============================================================================#
class Gait2392GeomCustomiser(object):
    gfield_disc = (8,8)
    convert_mm_to_m = True

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
        self.LL = None  # left lowerlimb object
        self.LR = None  # right lowerlimb object
        self.osimmodel = None  # opensim model

        self._init_osim_model()

    def _init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)

    def _check_geom_path(self):
        """
        Check that the directory for geom meshes exists. If not, create it.
        """
        geom_dir = os.path.join(self.config['osim_output_dir'], GEOM_DIR)
        if not os.path.isdir(geom_dir):
            os.mkdir(geom_dir)
        
    def set_left_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """

        self.LL = bonemodels.LowerLimbLeftAtlas('left lower limb')
        for gname, g in gfieldsdict.items():
            self.LL.set_bone_gfield(gname, g)
            self.LL.models[gname].update_acs()
        # self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
        # self.LL.set_bone_gfield('femur-left', gfieldsdict['femur-left'])
        # self.LL.set_bone_gfield('femur-left', gfieldsdict['femur-left'])
        # self.LL.set_bone_gfield('patella-left', gfieldsdict['patella-left'])
        # self.LL.set_bone_gfield('patella-right', gfieldsdict['patella-right'])
        # self.LL.set_bone_gfield('tibiafibula-left', gfieldsdict['tibiafibula-left'])
        # self.LL.set_bone_gfield('tibiafibula-right', gfieldsdict['tibiafibula-right'])

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
        if self.convert_mm_to_m:
            v_local *= 1e-3
        vtkwriter = vtktools.Writer(
                        v=v_local,
                        f=f,
                        filename=filename,
                        )
        vtkwriter.writeVTP()

    def cust_osim_pelvis(self):
        pelvis = self.LL.models['pelvis']
        osim_pelvis = self.osimmodel.bodies[OSIM_BODY_NAME_MAP['pelvis']]

        # update mass and inertial

        # update ground-pelvis joint
        pelvis_origin = pelvis.acs.o  
        self.osimmodel.joints['ground_pelvis'].locationInParent = pelvis_origin # in ground CS
        self.osimmodel.joints['ground_pelvis'].location =  (0,0,0) # in pelvis CS
        if self.convert_mm_to_m:
            self.osimmodel.joints['ground_pelvis'].locationInParent *= 1e-3  
            self.osimmodel.joints['ground_pelvis'].location *= 1e-3  

        ## pelvis_tilt

        ## pelvis_list

        ## pelvis_rotation

        # update mesh
        lhgf, sacgf, rhgf = _splitPelvisGFs(self.LL.models['pelvis'].gf)
        self._check_geom_path()

        ## sacrum.vtp
        sac_vtp_full_path = os.path.join(self.config['osim_output_dir'], GEOM_DIR, SACRUM_FILENAME)
        sac_vtp_osim_path = os.path.join(GEOM_DIR, SACRUM_FILENAME)
        self._save_vtp(sacgf, sac_vtp_full_path, pelvis.acs.map_local)

        ## pelvis.vtp
        rh_vtp_full_path = os.path.join(self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_RIGHT_FILENAME)
        rh_vtp_osim_path = os.path.join(GEOM_DIR, HEMIPELVIS_RIGHT_FILENAME)
        self._save_vtp(rhgf, rh_vtp_full_path, pelvis.acs.map_local)

        ## l_pelvis.vtp
        lh_vtp_full_path = os.path.join(self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_LEFT_FILENAME)
        lh_vtp_osim_path = os.path.join(GEOM_DIR, HEMIPELVIS_LEFT_FILENAME)
        self._save_vtp(lhgf, lh_vtp_full_path, pelvis.acs.map_local)

        osim_pelvis.setDisplayGeometryFileName(
            [sac_vtp_osim_path, rh_vtp_osim_path, lh_vtp_osim_path]
            )

    def cust_osim_femur_left(self):
        femur = self.LL.models['femur']
        pelvis = self.LL.models['pelvis']
        osim_femur = self.osimmodel.bodies[OSIM_BODY_NAME_MAP['femur-left']]

        # update hip_l joint
        lhjc = pelvis.landmarks['pelvis-LHJC']
        self.osimmodel.joints['hip_l'].locationInParent = pelvis.acs.map_local(lhjc[np.newaxis])[0]
        self.osimmodel.joints['hip_l'].location = femur.acs.map_local(lhjc[np.newaxis])[0]
        if self.convert_mm_to_m:
            self.osimmodel.joints['hip_l'].locationInParent *= 1e-3
            self.osimmodel.joints['hip_l'].location *= 1e-3

        # update coordinate defaults

        ## hip_flexion_l

        ## hip_adduction_l

        ## hip_rotation_l

        # update mesh l_femur.vtp
        self._check_geom_path()
        femur_vtp_full_path = os.path.join(self.config['osim_output_dir'], GEOM_DIR, FEMUR_LEFT_FILENAME)
        femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_LEFT_FILENAME)
        self._save_vtp(femur.gf, femur_vtp_full_path, femur.acs.map_local)
        osim_femur.setDisplayGeometryFileName([femur_vtp_osim_path,])

    def cust_osim_tibiafibula_left(self):
        tibfib = self.LL.models['tibiafibula']
        femur = self.LL.models['femur']
        osim_tibfib = self.osimmodel.bodies[OSIM_BODY_NAME_MAP['tibiafibula-left']]

        # update knee_l joint
        kjc = 0.5*(femur.landmarks['femur-MEC'] + femur.landmarks['femur-LEC'])
        self.osimmodel.joints['knee_l'].locationInParent = femur.acs.map_local(kjc[np.newaxis])[0]
        self.osimmodel.joints['knee_l'].location = tibfib.acs.map_local(kjc[np.newaxis])[0]
        if self.convert_mm_to_m:
            self.osimmodel.joints['knee_l'].locationInParent *= 1e-3
            self.osimmodel.joints['knee_l'].location *= 1e-3

         # update mesh
        tibgf, fibgf = _splitTibiaFibulaGFs(self.LL.models['tibiafibula'].gf)
        self._check_geom_path()

        # update mesh l_tibia.vtp
        self._check_geom_path()
        tib_vtp_full_path = os.path.join(
            self.config['osim_output_dir'],
            GEOM_DIR,
            TIBIA_LEFT_FILENAME
            )
        tib_vtp_osim_path = os.path.join(
            GEOM_DIR,
            TIBIA_LEFT_FILENAME
            )
        self._save_vtp(tibgf, tib_vtp_full_path, tibfib.acs.map_local)

        # update mesh l_fibula.vtp
        fib_vtp_full_path = os.path.join(
            self.config['osim_output_dir'],
            GEOM_DIR,
            FIBULA_LEFT_FILENAME
            )
        fib_vtp_osim_path = os.path.join(
            GEOM_DIR,
            FIBULA_LEFT_FILENAME
            )
        self._save_vtp(fibgf, fib_vtp_full_path, tibfib.acs.map_local)
        
        osim_tibfib.setDisplayGeometryFileName(
            [tib_vtp_osim_path, fib_vtp_osim_path]
            )

    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(self.config['osim_output_dir'], OSIM_FILENAME)
            )