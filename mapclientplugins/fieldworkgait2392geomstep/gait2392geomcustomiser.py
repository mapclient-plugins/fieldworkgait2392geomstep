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
from gias2.musculoskeletal.bonemodels import lowerlimbatlas
from gias2.musculoskeletal import osim
from gias2.musculoskeletal import fw_model_landmarks as fml

from transforms3d.euler import mat2euler

import opensim
import scaler

#=============================================================================#
SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
TEMPLATE_OSIM_PATH = os.path.join(SELF_DIR, 'data', 'gait2392_simbody.osim')
OSIM_FILENAME = 'gait2392_simbody.osim'
OSIM_BODY_NAME_MAP = {'pelvis': 'pelvis',
                      'femur-l': 'femur_l',
                      'femur-r': 'femur_r',
                      'tibiafibula-l': 'tibia_l',
                      'tibiafibula-r': 'tibia_r',
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
# SIDES = ('left', 'right', 'both')

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

def calc_pelvis_ground_angles(pelvis):
    """
    returns pelvis tilt, list, rotation relative to ground
    """
    globalCS = np.array(
        [[0,0,0],
         [0,0,1],
         [1,0,0],
         [0,1,0],
         ])
    pelvisACS = pelvis.acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(globalCS, pelvisACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    _list, tilt, rot = mat2euler(R, 'szxy')

    return -tilt, -_list, -rot

def calc_hip_angles(pelvis, femur, side):
    """
    returns hip flexion, adduction, rotation
    """
    pelvisACS = pelvis.acs.unit_array
    femurACS = femur.acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(pelvisACS, femurACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    rot, flex, add = mat2euler(R, 'szxy')

    if side=='l':
        return -flex, -rot,  add
    else:
        return -flex,  rot, -add

def calc_knee_angles(femur, tibfib, side):
    """
    returns knee flexion, adduction, rotation
    """
    femurACS = femur.acs.unit_array
    tibfibACS = tibfib.acs.unit_array
    # calc rotation matrix mapping pelvis ACS to femur ACS
    R = transform3D.directAffine(femurACS, tibfibACS)[:3,:3]

    # calculate euler angles from rotation matrix 
    rot, flex, add = mat2euler(R, 'szxy')

    if side=='l':
        return -flex, rot, -add
    else:
        return -flex, -rot, add

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
        _ll.update_tibiafibula([a,0,0])
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
    ankle_offset = np.array([0., -0.01, 0.])
    back_offset = np.array([0., 0.01, 0.])

    def __init__(self, config, gfieldsdict=None, ll=None):
        """
        Class for customising the OpenSim Gait2392 model's bodies and joints.
        Customisation is based on either an input LowerLimbAtlas instance or
        a dictionary of fieldwork geometric fields of each bone. Only one at
        most should be defined.

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
        gfieldsdict : dict [optional]
            Expected geometric field dict keys:
                pelvis
                femur-l
                femur-r
                patella-l
                patella-r
                tibiafibula-l
                tibiafibula-r
        ll : LowerLimbAtlas instance [optional]

        """
        self.config = config
        # self.ll_transform = None
        # self._pelvisRigid = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        # self._hipRot = np.array([0.0, 0.0, 0.0])
        # self._kneeRot = np.array([0.0, 0.0, 0.0])
        self.uniformScaling = 1.0
        self.pelvisScaling = 1.0
        self.femurScaling = 1.0
        self.petallaScaling = 1.0
        self.tibfibScaling = 1.0
        self.LL = None  # lowerlimb object
        self._hasInputLL = False 
        self.osimmodel = None  # opensim model
        self._unit_scaling = dim_unit_scaling(
                                self.config['in_unit'], self.config['out_unit']
                                )

        if gfieldsdict is not None:
            self.set_lowerlimb_gfields(gfieldsdict)
        if ll is not None:
            self.set_lowerlimb_atlas(ll)

    def init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)
        self._osimmodel_init_state = self.osimmodel._model.initSystem()

    def _check_geom_path(self):
        """
        Check that the directory for geom meshes exists. If not, create it.
        """
        geom_dir = os.path.join(self.config['osim_output_dir'], GEOM_DIR)
        if not os.path.isdir(geom_dir):
            os.mkdir(geom_dir)
    
    def set_lowerlimb_atlas(self, ll):
        self.LL = ll
        self._hasInputLL = True

        update_femur_opensim_acs(self.LL.ll_l.models['femur'])
        update_tibiafibula_opensim_acs(self.LL.ll_l.models['tibiafibula'])
        update_femur_opensim_acs(self.LL.ll_r.models['femur'])
        update_tibiafibula_opensim_acs(self.LL.ll_r.models['tibiafibula'])

    def set_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """
        self.set_2side_lowerlimb_gfields(gfieldsdict)

        # if self.config['side']=='left':
        #     self.set_left_lowerlimb_gfields(gfieldsdict)
        # elif self.config['side']=='right':
        #     self.set_right_lowerlimb_gfields(gfieldsdict)
        # elif self.config['side']=='both':
        #     self.set_2side_lowerlimb_gfields(gfieldsdict)

    # def set_left_lowerlimb_gfields(self, gfieldsdict):
    #     """
    #     Instantiate the lower limb object using input models
    #     """
    #     self.LL = bonemodels.LowerLimbLeftAtlas('left lower limb')
    #     self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
    #     self.LL.set_bone_gfield('femur', gfieldsdict['femur'])
    #     self.LL.set_bone_gfield('patella', gfieldsdict['patella'])
    #     self.LL.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula'])
    #     self.LL.models['pelvis'].update_acs()
    #     update_femur_opensim_acs(self.LL.models['femur'])
    #     update_tibiafibula_opensim_acs(self.LL.models['tibiafibula'])

    # def set_right_lowerlimb_gfields(self, gfieldsdict):
    #     """
    #     Instantiate the lower limb object using input models
    #     """
    #     self.LL = bonemodels.LowerLimbRightAtlas('right lower limb')
    #     self.LL.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
    #     self.LL.set_bone_gfield('femur', gfieldsdict['femur'])
    #     self.LL.set_bone_gfield('patella', gfieldsdict['patella'])
    #     self.LL.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula'])
    #     self.LL.models['pelvis'].update_acs()
    #     update_femur_opensim_acs(self.LL.models['femur'])
    #     update_tibiafibula_opensim_acs(self.LL.models['tibiafibula'])

    def set_2side_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """

        # left
        if not self._hasInputLL:
            ll_l = bonemodels.LowerLimbLeftAtlas('left lower limb')
            ll_l.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            ll_l.set_bone_gfield('femur', gfieldsdict['femur-l'])
            ll_l.set_bone_gfield('patella', gfieldsdict['patella-l'])
            ll_l.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-l'])
        else:
            ll_l = self.LL.ll_l
            if 'pelvis' in gfieldsdict:
                ll_l.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            if 'femur-l' in gfieldsdict:
                ll_l.set_bone_gfield('femur', gfieldsdict['femur-l'])
            if 'patella-l' in gfieldsdict:
                ll_l.set_bone_gfield('patella', gfieldsdict['patella-l'])
            if 'tibiafibula-l' in gfieldsdict:
                ll_l.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-l'])

        update_femur_opensim_acs(ll_l.models['femur'])
        update_tibiafibula_opensim_acs(ll_l.models['tibiafibula'])

        # right
        if not self._hasInputLL:
            ll_r = bonemodels.LowerLimbLeftAtlas('right lower limb')
            ll_r.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            ll_r.set_bone_gfield('femur', gfieldsdict['femur-r'])
            ll_r.set_bone_gfield('patella', gfieldsdict['patella-r'])
            ll_r.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-r'])
        else:
            ll_r = self.LL.ll_r
            if 'pelvis' in gfieldsdict:
                ll_r.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            if 'femur-r' in gfieldsdict:
                ll_r.set_bone_gfield('femur', gfieldsdict['femur-r'])
            if 'patella-r' in gfieldsdict:
                ll_r.set_bone_gfield('patella', gfieldsdict['patella-r'])
            if 'tibiafibula-r' in gfieldsdict:
                ll_r.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-r'])
        
        update_femur_opensim_acs(ll_r.models['femur'])
        update_tibiafibula_opensim_acs(ll_r.models['tibiafibula'])

        # 2side
        if not self._hasInputLL:
            self.LL = lowerlimbatlas.LowerLimbAtlas('lower limb')
            self.LL.ll_l = ll_l
            self.LL.ll_r = ll_r

        self.LL._update_model_dict()

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

        # update inertial and model muscle properties
        pelvis_sf = scaler.calc_pelvis_scale_factors(
                        self.LL, self._unit_scaling,
                    )
        pelvis_scaling = osim.Scale(
                            pelvis_sf,
                            'pelvis_scale',
                            'pelvis',
                        )
        self.osimmodel.scale(self._osimmodel_init_state, pelvis_scaling)

        # update ground-pelvis joint
        pelvis_origin = pelvis.acs.o  
        self.osimmodel.joints['ground_pelvis'].locationInParent = \
            pelvis_origin*self._unit_scaling # in ground CS
        self.osimmodel.joints['ground_pelvis'].location = \
            np.array((0,0,0), dtype=float)*self._unit_scaling  # in pelvis CS

        # update coordinate defaults
        pelvis_ground_joint = self.osimmodel.joints['ground_pelvis']
        if self._hasInputLL:
            tilt, _list, rot = self.LL.pelvis_rigid[3:]
        else:
            tilt, _list, rot = calc_pelvis_ground_angles(pelvis)

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

    def cust_osim_femur_l(self):
        self._cust_osim_femur('l')

    def cust_osim_femur_r(self):
        self._cust_osim_femur('r')

    def _cust_osim_femur(self, side):
        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        femur = self.LL.models['femur-'+side]
        pelvis = self.LL.models['pelvis']
        osim_femur = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP[
                            'femur-'+side
                            ]
                        ]

        # Apply general scaling to update inertial and muscle model properties
        femur_sf = scaler.calc_femur_scale_factors(
                    self.LL, self._unit_scaling, side
                    )
        femur_scaling = osim.Scale(
                            femur_sf,
                            'femur_{}_scale'.format(side),
                            'femur_{}'.format(side),
                        )
        self.osimmodel.scale(self._osimmodel_init_state, femur_scaling)

        # update hip joint
        if side=='l':
            hjc = pelvis.landmarks['pelvis-LHJC']
        else:
            hjc = pelvis.landmarks['pelvis-RHJC']
        self.osimmodel.joints['hip_{}'.format(side)].locationInParent = \
            pelvis.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling
        self.osimmodel.joints['hip_{}'.format(side)].location = \
            femur.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling

        # update coordinate defaults
        if self._hasInputLL:
            if side=='l':
                flex, rot, add = self.LL.hip_rot_l
            else:
                flex, rot, add = self.LL.hip_rot_r
            
        else:
            flex, rot, add = calc_hip_angles(pelvis, femur, side)

        hip_joint = self.osimmodel.joints['hip_{}'.format(side)]
        ## hip_flexion_l
        hip_joint.coordSets['hip_flexion_{}'.format(side)].defaultValue = flex
        ## hip_adduction_l
        hip_joint.coordSets['hip_adduction_{}'.format(side)].defaultValue = add
        ## hip_rotation_l
        hip_joint.coordSets['hip_rotation_{}'.format(side)].defaultValue = rot

        # update mesh l_femur.vtp
        self._check_geom_path()
        if side=='l':
            femur_vtp_full_path = os.path.join(
                self.config['osim_output_dir'], GEOM_DIR, FEMUR_LEFT_FILENAME
                )
            femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_LEFT_FILENAME)
        elif side=='r':
            femur_vtp_full_path = os.path.join(
                self.config['osim_output_dir'], GEOM_DIR, FEMUR_RIGHT_FILENAME
                )
            femur_vtp_osim_path = os.path.join(GEOM_DIR, FEMUR_RIGHT_FILENAME)

        self._save_vtp(femur.gf, femur_vtp_full_path, femur.acs.map_local)
        osim_femur.setDisplayGeometryFileName([femur_vtp_osim_path,])

    def _get_osim_knee_spline_xk(self, side):
        """
        Get the SimmSpline x values from the translation functions
        of the gati2392 knee
        """
        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        if side=='l':
            kj = self.osimmodel.joints['knee_l']
        else:
            kj = self.osimmodel.joints['knee_r']

        t1x = kj.getSimmSplineParams('translation1')[0]
        t2x = kj.getSimmSplineParams('translation2')[0]
        return t1x, t2x

    def _set_osim_knee_spline_xyk(self, x, y, side):
        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        if side=='l':
            kj = self.osimmodel.joints['knee_l']
        else:
            kj = self.osimmodel.joints['knee_r']

        kj.updateSimmSplineParams('translation1', x[0], y[0])
        kj.updateSimmSplineParams('translation2', x[1], y[1])

    def cust_osim_tibiafibula_l(self):
        self._cust_osim_tibiafibula('l')

    def cust_osim_tibiafibula_r(self):
        self._cust_osim_tibiafibula('r')

    def _cust_osim_tibiafibula(self, side):
        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        tibfib = self.LL.models['tibiafibula-'+side]
        femur = self.LL.models['femur-'+side]
        osim_tibfib = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP['tibiafibula-'+side]
                        ]

        # Apply general scaling to update inertial and muscle model properties
        tibia_sf = scaler.calc_tibia_scale_factors(
                    self.LL, self._unit_scaling, side
                    )
        tibia_scaling = osim.Scale(
                            tibia_sf,
                            'tibia_{}_scale'.format(side),
                            'tibia_{}'.format(side),
                        )
        self.osimmodel.scale(self._osimmodel_init_state, tibia_scaling)

        # update knee_l joint
        kjc = 0.5*(femur.landmarks['femur-MEC'] + femur.landmarks['femur-LEC'])
        tpc = 0.5*(tibfib.landmarks['tibiafibula-MC'] + tibfib.landmarks['tibiafibula-LC'])
        _d = -np.sqrt(((kjc - tpc)**2.0).sum())

        # Knee trans spline params are relative to the femoral head origin
        self.osimmodel.joints['knee_{}'.format(side)].locationInParent = \
            np.array([0,0,0], dtype=float)*self._unit_scaling 
        self.osimmodel.joints['knee_{}'.format(side)].location = \
            np.array([0,0,0], dtype=float)*self._unit_scaling

        # Knee spline values
        # get spline xk from osim files
        knee_spline_xk_1, knee_spline_xk_2 = self._get_osim_knee_spline_xk(side)
        knee_spline_xk = [knee_spline_xk_1, knee_spline_xk_2]
        # evaluate tib coord at xks
        if side=='l':
            knee_spline_yk_1 = _calc_knee_spline_coords(self.LL.ll_l, knee_spline_xk_1)*self._unit_scaling
            knee_spline_yk_2 = _calc_knee_spline_coords(self.LL.ll_l, knee_spline_xk_2)*self._unit_scaling
        else:
            knee_spline_yk_1 = _calc_knee_spline_coords(self.LL.ll_r, knee_spline_xk_1)*self._unit_scaling
            knee_spline_yk_2 = _calc_knee_spline_coords(self.LL.ll_r, knee_spline_xk_2)*self._unit_scaling
        knee_spline_yk = [knee_spline_yk_1[:,0], knee_spline_yk_2[:,1]]
        # set new spline yks
        self._set_osim_knee_spline_xyk(knee_spline_xk, knee_spline_yk, side)

        # Set input knee angle
        knee_joint = self.osimmodel.joints['knee_{}'.format(side)]
        if self._hasInputLL:
            if side=='l':
                flex, rot, add = self.LL._knee_rot_l
            else:
                flex, rot, add = self.LL._knee_rot_r
        else:
            flex, rot, add = calc_knee_angles(femur, tibfib, side)

        ## hip_flexion_l
        knee_joint.coordSets['knee_angle_{}'.format(side)].defaultValue = flex
        
         # update mesh
        tibgf, fibgf = _splitTibiaFibulaGFs(self.LL.models['tibiafibula-'+side].gf)
        self._check_geom_path()

        # update mesh l_tibia.vtp
        if side=='l':
            tibia_filename = TIBIA_LEFT_FILENAME
        if side=='r':
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
        if side=='l':
            fibula_filename = FIBULA_LEFT_FILENAME
        if side=='r':
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

    def cust_osim_ankle_l(self):
        # self._cust_osim_ankle('l')
        self._cust_osim_foot('l')

    def cust_osim_ankle_r(self):
        # self._cust_osim_ankle('r')
        self._cust_osim_foot('r')

    def _cust_osim_foot(self, side):
        """
        Customises foot models by applying opensim scaling to the foot segments,
        joints, and muscle sites.

        Segment topology in the foot is
        tibia -> ankle(j) -> talus -> subtalar(j) -> calcaneus -> mtp(j) -> toes
        """

        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        tibfib = self.LL.models['tibiafibula-'+side]
        femur = self.LL.models['femur-'+side]

        # Apply general scaling to update inertial and muscle model properties
        foot_sf = scaler.calc_whole_body_scale_factors(
                    self.LL, self._unit_scaling,
                    )
        print('foot {} scale factor: {}'.format(side, foot_sf))
        foot_scalings = [
            opysim.Scale(foot_sf, 'talus_{}_scaling'.format(side), 'talus_{}'.format(side)),
            opysim.Scale(foot_sf, 'calcn_{}_scaling'.format(side), 'calcn_{}'.format(side)),
            opysim.Scale(foot_sf, 'toes_{}_scaling'.format(side), 'toes_{}'.format(side)),
        ]
        self.osimmodel.scale(self._osimmodel_init_state, *foot_scalings)

        # set ankle joint parent location in custom tibiafibula
        ankle_centre = 0.5*(
            tibfib.landmarks['tibiafibula-MM'] + tibfib.landmarks['tibiafibula-LM']
            )
        self.osimmodel.joints['ankle_{}'.format(side)].locationInParent = \
            (tibfib.acs.map_local(ankle_centre[np.newaxis]).squeeze()*self._unit_scaling)+\
            self.ankle_offset

        # get scaling factor from femur size
        # scale_factor_0 = _calc_scaling_femur_length(femur)
        # scale_factor_array = np.array([scale_factor_0,]*3)  # isotropic scaling for now

        # scale_factor_array = scaler.calc_whole_body_scale_factors(self.LL, self._unit_scaling)
        
        # print('foot {} scale factor: {}'.format(side, scale_factor_array))

        # # scale bodies
        # _scale_body(self.osimmodel.bodies['talus_{}'.format(side)], scale_factor_array)
        # _scale_body(self.osimmodel.bodies['calcn_{}'.format(side)], scale_factor_array)
        # _scale_body(self.osimmodel.bodies['toes_{}'.format(side)], scale_factor_array)
        
        # # scale joints
        # ankle_scales = [
        #     osim.Scale(
        #         [scale_factor_0,]*3,
        #         'ankle_{}_scale_1'.format(side),
        #         'talus_{}'.format(side),
        #         ),
        #     ]
        # self.osimmodel.joints['ankle_{}'.format(side)].scale(*ankle_scales)

        # subtalar_scales = [
        #     osim.Scale(
        #         [scale_factor_0,]*3,
        #         'subtalar_{}_scale_0'.format(side),
        #         'talus_{}'.format(side),
        #         ),
        #     osim.Scale(
        #         [scale_factor_0,]*3,
        #         'subtalar_{}_scale_1'.format(side),
        #         'calcn_{}'.format(side),
        #         ),
        #     ]
        # self.osimmodel.joints['subtalar_{}'.format(side)].scale(*subtalar_scales)

        # mtp_scales = [
        #     osim.Scale(
        #         [scale_factor_0,]*3,
        #         'mtp_{}_scale_0'.format(side),
        #         'calcn_{}'.format(side),
        #         ),
        #     osim.Scale(
        #         [scale_factor_0,]*3,
        #         'mtp_{}_scale_1'.format(side),
        #         'toes_{}'.format(side),
        #         ),
        #     ]
        # self.osimmodel.joints['mtp_{}'.format(side)].scale(*mtp_scales)

        # # scale muscles of the foot
        # muscle_scales = [
        #     osim.Scale(scale_factor_array, 'mus_scale_0', 'talus_{}'.format(side)),
        #     osim.Scale(scale_factor_array, 'mus_scale_1', 'caln_{}'.format(side)),
        #     osim.Scale(scale_factor_array, 'mus_scale_2', 'toes_{}'.format(side)),
        #     ]

        # foot_muscles = _get_foot_muscles(self.osimmodel, side)

        # for mus in foot_muscles:
        #     print('scaling foot muscle: {}'.format(mus.name))
        #     mus.scale(self._osimmodel_init_state, *muscle_scales)

    def cust_torso(self):

        pelvis = self.LL.models['pelvis']
        femur_l = self.LL.models['femur-l']
        femur_r = self.LL.models['femur-r']

        # Apply general scaling to update inertial and muscle model properties
        torso_sf = scaler.calc_whole_body_scale_factors(
                    self.LL, self._unit_scaling,
                    )
        torso_scaling = osim.Scale(
                            torso_sf,
                            'torso_scale',
                            'torso',
                        )
        self.osimmodel.scale(self._osimmodel_init_state, torso_scaling)

        # set back joint parent location in custom pelvis
        sacrum_top = pelvis.landmarks['pelvis-SacPlat']
        self.osimmodel.joints['back'].locationInParent = \
            (pelvis.acs.map_local(sacrum_top[np.newaxis]).squeeze()*self._unit_scaling)+\
            self.back_offset

        # get scaling factor from femur size
        # scale_factor_fl = _calc_scaling_femur_length(femur_r)
        # scale_factor_fr = _calc_scaling_femur_length(femur_l)
        # scale_factor = 0.5*(scale_factor_fl + scale_factor_fr)
        # scale_factor_array = np.array([scale_factor,]*3)  # isotropic scaling for now

        # scale_factor_array = scaler.calc_whole_body_scale_factors(self.LL, self._unit_scaling)
        # print('torso scale factor: {}'.format(scale_factor_array))

        # # scale bodies
        # _scale_body(self.osimmodel.bodies['torso'], scale_factor_array)
        
        # # scale joints
        # back_scales = [
        #     osim.Scale(
        #         scale_factor_array,
        #         'back_scale_1',
        #         'torso',
        #         ),
        #     ]
        # self.osimmodel.joints['back'].scale(*back_scales)

        # # scale muscles of the torso
        # muscle_scales = [
        #     osim.Scale(scale_factor_array, 'mus_scale_0', 'torso'),
        #     ]
        # torso_muscles = _get_segment_muscles(self.osimmodel, 'torso')
        # for mus in torso_muscles:
        #     print('scaling torso muscle: {}'.format(mus.name))
        #     mus.scale(self._osimmodel_init_state, *muscle_scales)

    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(str(self.config['osim_output_dir']), OSIM_FILENAME)
            )

    def customise(self):
        self.cust_osim_pelvis()
        self.cust_osim_femur_l()
        self.cust_osim_femur_r()
        self.cust_osim_tibiafibula_l()
        self.cust_osim_tibiafibula_r()
        self.cust_osim_ankle_l()
        self.cust_osim_ankle_r()
        self.cust_torso()
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

def _get_foot_muscles(model, side):
    """
    Return osim muscles instances of muscles in the foot
    """
    foot_segs = set([
        'talus_{}'.format(side),
        'calcn_{}'.format(side),
        'toes_{}'.format(side),
        ])
    foot_muscles = []
    for mus in model.muscles.values():
        # check each path point to see if they are on a foot segment
        for pp in mus.getAllPathPoints():
            if pp.body.name in foot_segs:
                foot_muscles.append(mus)
                break

    return foot_muscles

def _get_segment_muscles(model, segname):
    """
    Return osim muscles instances of muscles in the defined segment
    """
    seg_muscles = []
    for mus in model.muscles.values():
        # check each path point to see if they are on the segment
        for pp in mus.getAllPathPoints():
            if pp.body.name==segname:
                seg_muscles.append(mus)
                break

    return seg_muscles

def _scale_body(body, sfarray):
    print('scaling {} by {}'.format(body.name, sfarray))
    body.scale(sfarray, False)
    body.scaleInertialProperties(sfarray, True)
    # body.scaleMass(sfarray.prod()) # mass should be cubed of scale factor