"""
Gait2392 customisation
"""
import os
import numpy as np
import copy

from gias2.common import transform3D
from gias2.fieldwork.field import geometric_field
from gias2.mesh import vtktools
from gias2.musculoskeletal import mocap_landmark_preprocess
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfit
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfitscaling
from gias2.musculoskeletal import osim

from transforms3d.euler import mat2euler

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

def update_femur_opensim_acs(femur_model):
    femur_model.acs.update(
        *bonemodels.model_alignment.createFemurACSOpenSim(
            femur_model.landmarks['femur-HC'],
            femur_model.landmarks['femur-MEC'],
            femur_model.landmarks['femur-LEC'],
            side=femur_model.side
            )
        )
# bonemodels.FemurModel.update_acs = update_femur_opensim_acs

def update_tibiafibula_opensim_acs(tibiafibula_model):
    tibiafibula_model.acs.update(
        *bonemodels.model_alignment.createTibiaFibulaACSOpenSim(
            tibiafibula_model.landmarks['tibiafibula-MM'],
            tibiafibula_model.landmarks['tibiafibula-LM'],
            tibiafibula_model.landmarks['tibiafibula-MC'],
            tibiafibula_model.landmarks['tibiafibula-LC'],
            side=tibiafibula_model.side
            )
        )
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

def calc_pelvis_ground_angles(ll):
    """
    returns pelvis tilt, list, rotation relative to ground
    """
    globalCS = np.array(
        [[0,0,0],
         [0,0,1],
         [1,0,0],
         [0,1,0],
         ])
    pelvisACS = ll.models['pelvis'].acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(globalCS, pelvisACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    _list, tilt, rot = mat2euler(R, 'szxy')

    return -tilt, -_list, -rot

def calc_hip_angles(ll):
    """
    returns hip flexion, adduction, rotation
    """
    pelvisACS = ll.models['pelvis'].acs.unit_array
    femurACS = ll.models['femur'].acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(pelvisACS, femurACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    add, flex, rot = mat2euler(R, 'szxy')

    return -flex, -rot, -add

def calc_knee_angles(ll):
    """
    returns knee flexion, adduction, rotation
    """
    femurACS = ll.models['femur'].acs.unit_array
    tibfibACS = ll.models['tibiafibula'].acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(femurACS, tibfibACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    rot, flex, add = mat2euler(R, 'szxy')

    return -flex, -rot, -add

#=============================================================================#
class Gait2392GeomCustomiser(object):
    gfield_disc = (6,6)

    def __init__(self, config):
        self.config = config
        self.ll_transform = None
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
        self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
        self.LL.set_bone_gfield('femur', gfieldsdict['femur'])
        self.LL.set_bone_gfield('patella', gfieldsdict['patella'])
        self.LL.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula'])

        # for gname, g in gfieldsdict.items():
        #     self.LL.set_bone_gfield(gname, g)
        #     self.LL.models[gname].update_acs()

        self.LL.models['pelvis'].update_acs()
        update_femur_opensim_acs(self.LL.models['femur'])
        update_tibiafibula_opensim_acs(self.LL.models['tibiafibula'])

        # self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
        # self.LL.set_bone_gfield('femur-left', gfieldsdict['femur-left'])
        # self.LL.set_bone_gfield('femur-left', gfieldsdict['femur-left'])
        # self.LL.set_bone_gfield('patella-left', gfieldsdict['patella-left'])
        # self.LL.set_bone_gfield('patella-right', gfieldsdict['patella-right'])
        # self.LL.set_bone_gfield('tibiafibula-left', gfieldsdict['tibiafibula-left'])
        # self.LL.set_bone_gfield('tibiafibula-right', gfieldsdict['tibiafibula-right'])

    def _save_vtp(self, gf, filename, bodycoordmapper):
        v, f = gf.triangulate(self.gfield_disc)
        # f = f[:,::-1]
        v_local = bodycoordmapper(v)
        if self.config['convert_mm_to_m']:
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
        if self.config['convert_mm_to_m']:
            self.osimmodel.joints['ground_pelvis'].locationInParent *= 1e-3  
            self.osimmodel.joints['ground_pelvis'].location *= 1e-3  

        # update coordinate defaults
        pelvis_ground_joint = self.osimmodel.joints['ground_pelvis']
        if self.ll_transform is None:
            tilt, _list, rot = calc_pelvis_ground_angles(self.LL)
        else:
            tilt, _list, rot = self.ll_transform.pelvisRigid[3:]
        ## tilt
        pelvis_ground_joint.coordSets['pelvis_tilt'].defaultValue = tilt
        ## list
        pelvis_ground_joint.coordSets['pelvis_list'].defaultValue = _list
        ## rotation
        pelvis_ground_joint.coordSets['pelvis_rotation'].defaultValue = rot

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
        if self.config['convert_mm_to_m']:
            self.osimmodel.joints['hip_l'].locationInParent *= 1e-3
            self.osimmodel.joints['hip_l'].location *= 1e-3

        # update coordinate defaults
        hip_joint = self.osimmodel.joints['hip_l']
        if self.ll_transform is None:
            flex, rot, add = calc_hip_angles(self.LL)
        else:
            flex, rot, add = -1.0*self.ll_transform.hipRot
        ## hip_flexion_l
        hip_joint.coordSets['hip_flexion_l'].defaultValue = flex
        ## hip_adduction_l
        hip_joint.coordSets['hip_adduction_l'].defaultValue = add
        ## hip_rotation_l
        hip_joint.coordSets['hip_rotation_l'].defaultValue = rot

        # update mesh l_femur.vtp
        self._check_geom_path()
        femur_vtp_full_path = os.path.join(
            self.config['osim_output_dir'],
            GEOM_DIR, FEMUR_LEFT_FILENAME
            )
        femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_LEFT_FILENAME)
        self._save_vtp(femur.gf, femur_vtp_full_path, femur.acs.map_local)
        osim_femur.setDisplayGeometryFileName([femur_vtp_osim_path,])

    def cust_osim_tibiafibula_left(self):
        tibfib = self.LL.models['tibiafibula']
        femur = self.LL.models['femur']
        osim_tibfib = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP['tibiafibula-left']
                        ]

        # update knee_l joint
        kjc = 0.5*(femur.landmarks['femur-MEC'] + femur.landmarks['femur-LEC'])
        tpc = 0.5*(tibfib.landmarks['tibiafibula-MC'] + tibfib.landmarks['tibiafibula-LC'])
        _d = -np.sqrt(((kjc - tpc)**2.0).sum())
        # self.osimmodel.joints['knee_l'].locationInParent = femur.acs.map_local(
        #     kjc[np.newaxis]
        #     )[0] #- [0, -3.74557519e+02, 0]

        # not sure why, the femur origin is in the head, so the knee centre should not be 0,0,0
        # self.osimmodel.joints['knee_l'].locationInParent = [0,0,0] 
        kjc_femur = femur.acs.map_local(kjc[np.newaxis]).squeeze()
        kjc_femur[1] = -15 #_d
        self.osimmodel.joints['knee_l'].locationInParent = kjc_femur
        self.osimmodel.joints['knee_l'].location = tibfib.acs.map_local(
                                                    kjc[np.newaxis]
                                                    ).squeeze()
        if self.config['convert_mm_to_m']:
            self.osimmodel.joints['knee_l'].locationInParent *= 1e-3
            self.osimmodel.joints['knee_l'].location *= 1e-3

        # update coordinate defaults
        # SKIP FOR NOW SINCE WE NEED SPLINE PARAMETERS TOO
        
        knee_joint = self.osimmodel.joints['knee_l']
        if self.ll_transform is None:
            flex, rot, add = calc_knee_angles(self.LL)
        else:
            flex, rot, add = -1.0*self.ll_transform._kneeRot
        ## hip_flexion_l
        knee_joint.coordSets['knee_angle_l'].defaultValue = flex
        

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

    def cust_osim_ankle_left(self):
        tibfib = self.LL.models['tibiafibula']

        # very rough/bad estimate
        ankle_centre = 0.5*(tibfib.landmarks['tibiafibula-MM'] + tibfib.landmarks['tibiafibula-LM'])
        self.osimmodel.joints['ankle_l'].locationInParent = tibfib.acs.map_local(
                                                                ankle_centre[np.newaxis]
                                                                ).squeeze()
        self.osimmodel.joints['ankle_l'].location = [0,10,0]
        if self.config['convert_mm_to_m']:
            self.osimmodel.joints['ankle_l'].locationInParent *= 1e-3
            self.osimmodel.joints['ankle_l'].location *= 1e-3


    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(str(self.config['osim_output_dir']), OSIM_FILENAME)
            )

    def customise(self):
        self.cust_osim_pelvis()
        self.cust_osim_femur_left()
        self.cust_osim_tibiafibula_left()
        self.cust_osim_ankle_left()
        if self.config['write_osim_file']:
            self.write_cust_osim_model()