"""
MAP Client, a program to generate detailed musculoskeletal models for OpenSim.
    Copyright (C) 2012  University of Auckland
    
This file is part of MAP Client. (http://launchpad.net/mapclient)

    MAP Client is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MAP Client is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAP Client.  If not, see <http://www.gnu.org/licenses/>..
"""

"""
OpenSim Gait2392 bodies customisation
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
TEMPLATE_OSIM_PATH = os.path.join(SELF_DIR, 'data', 'gait2392_simbody.osim')
OSIM_FILENAME = 'gait2392_simbody.osim'
OSIM_BODY_NAME_MAP = {'pelvis': 'pelvis',
                      'femur-left': 'femur_l',
                      'femur-right': 'femur_r',
                      'tibiafibula-left': 'tibia_l',
                      'tibiafibula-right': 'tibia_r',
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

GEOM_DIR = 'geom'
SACRUM_FILENAME = 'sacrum.vtp'
HEMIPELVIS_RIGHT_FILENAME = 'pelvis.vtp'
HEMIPELVIS_LEFT_FILENAME = 'l_pelvis.vtp'
FEMUR_LEFT_FILENAME = 'l_femur.vtp'
TIBIA_LEFT_FILENAME = 'l_tibia.vtp'
FIBULA_LEFT_FILENAME = 'l_fibula.vtp'
FEMUR_RIGHT_FILENAME = 'r_femur.vtp'
TIBIA_RIGHT_FILENAME = 'r_tibia.vtp'
FIBULA_RIGHT_FILENAME = 'r_fibula.vtp'

VALID_UNITS = ('nm', 'um', 'mm', 'cm', 'm', 'km')
SIDES = ('left', 'right')

#=============================================================================#
def dim_unit_scaling(in_unit, out_unit):
    """
    Calculate the scaling factor to convert from the input unit (in_unit) to
    the output unit (out_unit). in_unit and out_unit must be a string and one
    of ['nm', 'um', 'mm', 'cm', 'm', 'km']. 

    inputs
    ======
    in_unit : str
        Input unit
    out_unit :str
        Output unit

    returns
    =======
    scaling_factor : float
    """

    unit_vals = {
        'nm': 1e-9,
        'um': 1e-6,
        'mm': 1e-3,
        'cm': 1e-2,
        'm':  1.0,
        'km': 1e3,
        }

    if in_unit not in unit_vals:
        raise ValueError(
            'Invalid input unit {}. Must be one of {}'.format(
                in_unit, list(unit_vals.keys())
                )
            )
    if out_unit not in unit_vals:
        raise ValueError(
            'Invalid input unit {}. Must be one of {}'.format(
                in_unit, list(unit_vals.keys())
                )
            )

    return unit_vals[in_unit]/unit_vals[out_unit]

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

def _calc_knee_spline_coords(ll, flex_angles):
    """
    Calculates the cubic spline values for the knee joint through specified
    angles. The values are the coordinates of the tibia frame origin relative
    to the femur frame.

    inputs
    ======
    ll : LowerLimbLeftAtlas instance
    flex_angles : 1d ndarray
        a list of n knee angles at which to sample tibia location relative to the
        femur ACS. Only flexion supported in 2392.

    returns
    =======
    y : n x 3 ndarray
        Array of the tibia frame origin relative to the femur frame at each
        knee angle.
    """

    _ll = copy.deepcopy(ll)
    # restore original ACSs
    _ll.models['femur'].update_acs()
    _ll.models['tibiafibula'].update_acs()
    # sample tibia ACS origin at each flexion angle
    tib_os = []
    for a in flex_angles:
        _ll.update_tibiafibula([a,])
        tib_o = 0.5*(_ll.models['tibiafibula'].landmarks['tibiafibula-LC'] + 
                     _ll.models['tibiafibula'].landmarks['tibiafibula-MC']
                     )
        tib_os.append(tib_o)

    update_femur_opensim_acs(_ll.models['femur'])
    y = _ll.models['femur'].acs.map_local(np.array(tib_os))
    # y = np.array([y[:,2], y[:,1], y[:,0]]).T # reverse dims
    return y

#=============================================================================#
class Gait2392GeomCustomiser(object):
    gfield_disc = (6,6)

    def __init__(self, config):
        """
        Class for customising the OpenSim Gait2392 model's bodies and joints.

        inputs
        ======
        config : dict
            Dict of configurable options:
            'osim_output_dir' : str
                Path to write out the customised .osim file.
            'write_osim_file' : bool
                If True, write customised .osim file to osim_output_dir.
            'in_unit' : str
                Input model's coordinate units
            'out_unit' : str
                Output model's coordinate units
            'side' : str
                Which limb to customised. Currently 'left' or 'right'.

        """
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
        self.LL = None  # lowerlimb object
        self.osimmodel = None  # opensim model
        self._unit_scaling = dim_unit_scaling(
                                self.config['in_unit'], self.config['out_unit']
                                )

        # self._init_osim_model()

    def init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)

    def _check_geom_path(self):
        """
        Check that the directory for geom meshes exists. If not, create it.
        """
        geom_dir = os.path.join(self.config['osim_output_dir'], GEOM_DIR)
        if not os.path.isdir(geom_dir):
            os.mkdir(geom_dir)
        
    def set_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """
        if self.config['side']=='left':
            self.set_left_lowerlimb_gfields(gfieldsdict)
        elif self.config['side']=='right':
            self.set_right_lowerlimb_gfields(gfieldsdict)

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

    def set_right_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """
        self.LL = bonemodels.LowerLimbRightAtlas('right lower limb')
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
        v_local *= self._unit_scaling
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
        self.osimmodel.joints['ground_pelvis'].locationInParent = \
            pelvis_origin*self._unit_scaling # in ground CS
        self.osimmodel.joints['ground_pelvis'].location = \
            np.array((0,0,0), dtype=float)*self._unit_scaling  # in pelvis CS

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
        sac_vtp_full_path = os.path.join(
            self.config['osim_output_dir'], GEOM_DIR, SACRUM_FILENAME
            )
        sac_vtp_osim_path = os.path.join(GEOM_DIR, SACRUM_FILENAME)
        self._save_vtp(sacgf, sac_vtp_full_path, pelvis.acs.map_local)

        ## pelvis.vtp
        rh_vtp_full_path = os.path.join(
            self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_RIGHT_FILENAME
            )
        rh_vtp_osim_path = os.path.join(GEOM_DIR, HEMIPELVIS_RIGHT_FILENAME)
        self._save_vtp(rhgf, rh_vtp_full_path, pelvis.acs.map_local)

        ## l_pelvis.vtp
        lh_vtp_full_path = os.path.join(
            self.config['osim_output_dir'], GEOM_DIR, HEMIPELVIS_LEFT_FILENAME
            )
        lh_vtp_osim_path = os.path.join(GEOM_DIR, HEMIPELVIS_LEFT_FILENAME)
        self._save_vtp(lhgf, lh_vtp_full_path, pelvis.acs.map_local)

        osim_pelvis.setDisplayGeometryFileName(
            [sac_vtp_osim_path, rh_vtp_osim_path, lh_vtp_osim_path]
            )

    def cust_osim_femur(self):
        femur = self.LL.models['femur']
        pelvis = self.LL.models['pelvis']
        osim_femur = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP[
                            'femur-{}'.format(
                                self.config['side']
                                )
                            ]
                        ]

        side_suffix = self.config['side'][0].lower()

        # update hip joint
        if self.config['side']=='left':
            hjc = pelvis.landmarks['pelvis-LHJC']
        elif self.config['side']=='right':
            hjc = pelvis.landmarks['pelvis-RHJC']
        self.osimmodel.joints['hip_{}'.format(side_suffix)].locationInParent = \
            pelvis.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling
        self.osimmodel.joints['hip_{}'.format(side_suffix)].location = \
            femur.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling

        # update coordinate defaults
        hip_joint = self.osimmodel.joints['hip_{}'.format(side_suffix)]
        if self.ll_transform is None:
            flex, rot, add = calc_hip_angles(self.LL)
        else:
            flex, rot, add = -1.0*self.ll_transform.hipRot
        ## hip_flexion_l
        hip_joint.coordSets['hip_flexion_{}'.format(side_suffix)].defaultValue = flex
        ## hip_adduction_l
        hip_joint.coordSets['hip_adduction_{}'.format(side_suffix)].defaultValue = add
        ## hip_rotation_l
        hip_joint.coordSets['hip_rotation_{}'.format(side_suffix)].defaultValue = rot

        # update mesh l_femur.vtp
        self._check_geom_path()
        if self.config['side']=='left':
            femur_vtp_full_path = os.path.join(
                self.config['osim_output_dir'], GEOM_DIR, FEMUR_LEFT_FILENAME
                )
            femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_LEFT_FILENAME)
        elif self.config['side']=='right':
            femur_vtp_full_path = os.path.join(
                self.config['osim_output_dir'], GEOM_DIR, FEMUR_RIGHT_FILENAME
                )
            femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_RIGHT_FILENAME)

        self._save_vtp(femur.gf, femur_vtp_full_path, femur.acs.map_local)
        osim_femur.setDisplayGeometryFileName([femur_vtp_osim_path,])

    def _get_osim_knee_spline_xk(self):
        """
        Get the SimmSpline x values from the translation functions
        of the gati2392 knee
        """
        if self.config['side']=='left':
            kj = self.osimmodel.joints['knee_l']
        elif self.config['side']=='right':
            kj = self.osimmodel.joints['knee_r']

        t1x = kj.getSimmSplineParams('translation1')[0]
        t2x = kj.getSimmSplineParams('translation2')[0]
        return t1x, t2x

    def _set_osim_knee_spline_xyk(self, x, y):
        if self.config['side']=='left':
            kj = self.osimmodel.joints['knee_l']
        elif self.config['side']=='right':
            kj = self.osimmodel.joints['knee_r']

        kj.updateSimmSplineParams('translation1', x[0], y[0])
        kj.updateSimmSplineParams('translation2', x[1], y[1])

    def cust_osim_tibiafibula(self):
        side_suffix = self.config['side'][0].lower()
        tibfib = self.LL.models['tibiafibula']
        femur = self.LL.models['femur']
        osim_tibfib = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP['tibiafibula-{}'.format(self.config['side'])]
                        ]

        # update knee_l joint
        kjc = 0.5*(femur.landmarks['femur-MEC'] + femur.landmarks['femur-LEC'])
        tpc = 0.5*(tibfib.landmarks['tibiafibula-MC'] + tibfib.landmarks['tibiafibula-LC'])
        _d = -np.sqrt(((kjc - tpc)**2.0).sum())

        # Knee trans spline params are relative to the femoral head origin
        self.osimmodel.joints['knee_{}'.format(side_suffix)].locationInParent = \
            np.array([0,0,0], dtype=float)*self._unit_scaling 
        self.osimmodel.joints['knee_{}'.format(side_suffix)].location = \
            np.array([0,0,0], dtype=float)*self._unit_scaling

        # Knee spline values
        # get spline xk from osim files
        knee_spline_xk_1, knee_spline_xk_2 = self._get_osim_knee_spline_xk()
        knee_spline_xk = [knee_spline_xk_1, knee_spline_xk_2]
        # evaluate tib coord at xks
        knee_spline_yk_1 = _calc_knee_spline_coords(self.LL, knee_spline_xk_1)*self._unit_scaling
        knee_spline_yk_2 = _calc_knee_spline_coords(self.LL, knee_spline_xk_2)*self._unit_scaling
        knee_spline_yk = [knee_spline_yk_1[:,0], knee_spline_yk_2[:,1]]
        # set new spline yks
        self._set_osim_knee_spline_xyk(knee_spline_xk, knee_spline_yk)

        # Set input knee angle
        knee_joint = self.osimmodel.joints['knee_{}'.format(side_suffix)]
        if self.ll_transform is None:
            flex, rot, add = calc_knee_angles(self.LL)
        else:
            flex, rot, add = self.ll_transform._kneeRot
        ## hip_flexion_l
        knee_joint.coordSets['knee_angle_{}'.format(side_suffix)].defaultValue = flex
        
         # update mesh
        tibgf, fibgf = _splitTibiaFibulaGFs(self.LL.models['tibiafibula'].gf)
        self._check_geom_path()

        # update mesh l_tibia.vtp
        if side_suffix=='l':
            tibia_filename = TIBIA_LEFT_FILENAME
        if side_suffix=='r':
            tibia_filename = TIBIA_RIGHT_FILENAME
        self._check_geom_path()
        tib_vtp_full_path = os.path.join(
            self.config['osim_output_dir'],
            GEOM_DIR,
            tibia_filename,
            )
        tib_vtp_osim_path = os.path.join(
            GEOM_DIR,
            tibia_filename,
            )
        self._save_vtp(tibgf, tib_vtp_full_path, tibfib.acs.map_local)

        # update mesh l_fibula.vtp
        if side_suffix=='l':
            fibula_filename = FIBULA_LEFT_FILENAME
        if side_suffix=='r':
            fibula_filename = FIBULA_RIGHT_FILENAME
        fib_vtp_full_path = os.path.join(
            self.config['osim_output_dir'],
            GEOM_DIR,
            fibula_filename,
            )
        fib_vtp_osim_path = os.path.join(
            GEOM_DIR,
            fibula_filename,
            )
        self._save_vtp(fibgf, fib_vtp_full_path, tibfib.acs.map_local)
        
        osim_tibfib.setDisplayGeometryFileName(
            [tib_vtp_osim_path, fib_vtp_osim_path]
            )

    def cust_osim_ankle(self):
        side_suffix = self.config['side'][0].lower()
        tibfib = self.LL.models['tibiafibula']

        # very rough/bad estimate
        ankle_centre = 0.5*(tibfib.landmarks['tibiafibula-MM'] + tibfib.landmarks['tibiafibula-LM'])
        self.osimmodel.joints['ankle_{}'.format(side_suffix)].locationInParent = \
            tibfib.acs.map_local(ankle_centre[np.newaxis]).squeeze()*self._unit_scaling
        self.osimmodel.joints['ankle_{}'.format(side_suffix)].location = \
            np.array([0,10,0], dtype=float)*self._unit_scaling

    def cust_gen_scaling(self):
        """
        Calculate a general isotropic scaling factor to apply to non-patient
        specific segments
        """

        # calculate scaling factor from femoral length.
        # f = new_femur_length / opensim_femur_length
        scaling_factor = _calc_scaling_femur_length(self.LL.models['femur'])

        # apply scaling factor
        _apply_uniform_scaling(self.osimmodel.bodies['talus_r'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['calcn_r'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['toes_r'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['talus_l'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['calcn_l'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['toes_l'], scaling_factor)
        _apply_uniform_scaling(self.osimmodel.bodies['torso'], scaling_factor)


    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(str(self.config['osim_output_dir']), OSIM_FILENAME)
            )

    def customise(self):
        self.cust_osim_pelvis()
        self.cust_osim_femur()
        self.cust_osim_tibiafibula()
        self.cust_osim_ankle()
        if self.config['write_osim_file']:
            self.write_cust_osim_model()

def _calc_scaling_femur_length(new_femur):
    """
    Calculate the ratio l_new / l_osim where l is the distance between
    a femur's femoral head and mid point of the epicondyles.

    for the gait2392 femur:
    MEC = [ 0.01173 , -0.3978  ,  0.036877]
    LEC = [-0.002143, -0.406302, -0.037033]
    HC = [ 0.0004065,  0.0030935,  0.003309 ]
    """

    L_OSIM = 405.1824

    new_HC = fml.makeLandmarkEvaluator(
        'femur-HC', new_femur.gf
        )(new_femur.gf.field_parameters)
    new_MEC = fml.makeLandmarkEvaluator(
        'femur-MEC', new_femur.gf
        )(new_femur.gf.field_parameters)
    new_LEC = fml.makeLandmarkEvaluator(
        'femur-LEC', new_femur.gf
        )(new_femur.gf.field_parameters)
    
    l_new = np.sqrt(((new_HC - 0.5*(new_MEC + new_LEC))**2.0).sum())

    return l_new/L_OSIM