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

VALID_MODEL_MARKERS = sorted(list(scaler.virtualmarker.markers.keys()))

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
    # back_offset = np.array([0., 0.0, 0.])

    _body_scalers = {
                'torso': scaler.calc_whole_body_scale_factors,
                'pelvis': scaler.calc_pelvis_scale_factors,
                'femur_l': scaler.calc_femur_scale_factors,
                'femur_r': scaler.calc_femur_scale_factors,
                'tibia_l': scaler.calc_tibia_scale_factors,
                'tibia_r': scaler.calc_tibia_scale_factors,
                'talus_l': scaler.calc_whole_body_scale_factors,
                'talus_r': scaler.calc_whole_body_scale_factors,
                'calcn_l': scaler.calc_whole_body_scale_factors,
                'calcn_r': scaler.calc_whole_body_scale_factors,
                'toes_l': scaler.calc_whole_body_scale_factors,
                'toes_r': scaler.calc_whole_body_scale_factors,
                }

    def __init__(self, config, gfieldsdict=None, ll=None, verbose=True):
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
        self.uniform_scaling = 1.0
        self.pelvis_scaling = 1.0
        self.femur_scaling = 1.0
        self.petalla_scaling = 1.0
        self.tibfib_scaling = 1.0
        self.LL = None  # lowerlimb object
        self._hasInputLL = False 
        self.osimmodel = None  # opensim model
        self.markerset = None  # markerset associated with opensim model
        self.input_markers = {} # input marker name : input marker coords
        self.verbose = verbose
        self._unit_scaling = dim_unit_scaling(
                                self.config['in_unit'], self.config['out_unit']
                                )

        if gfieldsdict is not None:
            self.set_lowerlimb_gfields(gfieldsdict)
        if ll is not None:
            self.set_lowerlimb_atlas(ll)

        self._body_scale_factors = {}

    def init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)
        self._osimmodel_init_state = self.osimmodel._model.initSystem()
        self._original_segment_masses = dict([(b.name, b.mass) for b in self.osimmodel.bodies.values()])

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

    def _get_osimbody_scale_factors(self, bodyname):
        """
        Returns the scale factor for a body. Caches scale factors
        that have already been calculated

        inputs
        ------
        bodyname : str
            Gait2392 name of a body

        returns
        -------
        sf : length 3 ndarray
            scale factor array
        """
        
        if bodyname not in self._body_scale_factors:
            sf = self._body_scalers[bodyname](self.LL, self._unit_scaling)
            self._body_scale_factors[bodyname] = sf

        return self._body_scale_factors[bodyname]

    def cust_osim_pelvis(self):

        if self.verbose:
            print('\nCUSTOMISING PELVIS...')

        pelvis = self.LL.models['pelvis']
        osim_pelvis = self.osimmodel.bodies[OSIM_BODY_NAME_MAP['pelvis']]

        # scale inertial properties
        # sf = scaler.calc_pelvis_scale_factors(
        #         self.LL, self._unit_scaling,
        #         )
        sf = self._get_osimbody_scale_factors('pelvis')
        scaler.scale_body_mass_inertia(osim_pelvis, sf)
        if self.verbose:
            print('scale factor: {}'.format(sf))

        # update ground-pelvis joint
        if self.verbose:
            print('updating pelvis-ground joint...')

        pelvis_origin = pelvis.acs.o  
        self.osimmodel.joints['ground_pelvis'].locationInParent = \
            pelvis_origin*self._unit_scaling # in ground CS
        self.osimmodel.joints['ground_pelvis'].location = \
            np.array((0,0,0), dtype=float)*self._unit_scaling  # in pelvis CS

        if self.verbose:
            print(
                'location in parent: {}'.format(
                    self.osimmodel.joints['ground_pelvis'].locationInParent
                    )
                )
            print(
                'location: {}'.format(
                    self.osimmodel.joints['ground_pelvis'].location
                    )
                )

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

        if self.verbose:
            print(
                'pelvis tilt, list, rotation: {:5.2f}, {:5.2f}, {:5.2f}'.format(
                    pelvis_ground_joint.coordSets['pelvis_tilt'].defaultValue,
                    pelvis_ground_joint.coordSets['pelvis_list'].defaultValue,
                    pelvis_ground_joint.coordSets['pelvis_rotation'].defaultValue,
                    )
                )

        # update mesh
        if self.verbose:
            print('updating visual geometry...')

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

        if self.verbose:
            print('\nCUSTOMISING FEMUR {}'.format(side.upper()))

        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        femur = self.LL.models['femur-'+side]
        pelvis = self.LL.models['pelvis']
        osim_femur = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP[
                            'femur-'+side
                            ]
                        ]

        # scale inertial properties
        # sf = scaler.calc_femur_scale_factors(
        #         self.LL, self._unit_scaling,
        #         side=None,
        #         )
        sf = self._get_osimbody_scale_factors('femur_'+side)
        scaler.scale_body_mass_inertia(osim_femur, sf)
        if self.verbose:
            print('scale factor: {}'.format(sf))

        # remove multiplier functions from hip joint translations
        hip = self.osimmodel.joints['hip_{}'.format(side)]
        _remove_multiplier(hip.spatialTransform.get_translation1())
        _remove_multiplier(hip.spatialTransform.get_translation2())
        _remove_multiplier(hip.spatialTransform.get_translation3())

        # update hip joint
        if self.verbose:
            print('updating hip {} joint...'.format(side))

        if side=='l':
            hjc = pelvis.landmarks['pelvis-LHJC']
        else:
            hjc = pelvis.landmarks['pelvis-RHJC']
        self.osimmodel.joints['hip_{}'.format(side)].locationInParent = \
            pelvis.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling
        self.osimmodel.joints['hip_{}'.format(side)].location = \
            femur.acs.map_local(hjc[np.newaxis])[0] * self._unit_scaling

        if self.verbose:
            print(
                'location in parent: {}'.format(
                    self.osimmodel.joints['hip_{}'.format(side)].locationInParent
                    )
                )
            print(
                'location: {}'.format(
                    self.osimmodel.joints['hip_{}'.format(side)].location
                    )
                )

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

        if self.verbose:
            print(
                'hip flexion, adduction, rotation: {:5.2f}, {:5.2f}, {:5.2f}'.format(
                    hip_joint.coordSets['hip_flexion_{}'.format(side)].defaultValue,
                    hip_joint.coordSets['hip_adduction_{}'.format(side)].defaultValue,
                    hip_joint.coordSets['hip_rotation_{}'.format(side)].defaultValue,
                    )
                )

        # update mesh l_femur.vtp
        if self.verbose:
            print('updating visual geometry...')

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
        if self.verbose:
            print('\nCUSTOMISING TIBIA {}'.format(side.upper()))

        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        tibfib = self.LL.models['tibiafibula-'+side]
        femur = self.LL.models['femur-'+side]
        osim_tibfib = self.osimmodel.bodies[
                        OSIM_BODY_NAME_MAP['tibiafibula-'+side]
                        ]

        # scale inertial properties
        # sf = scaler.calc_tibia_scale_factors(
        #         self.LL, self._unit_scaling,
        #         side=None,
        #         )
        sf = self._get_osimbody_scale_factors('tibia_'+side)
        scaler.scale_body_mass_inertia(osim_tibfib, sf)
        if self.verbose:
            print('scale factor: {}'.format(sf))

        # recover knee joint simmspline
        knee = self.osimmodel.joints['knee_{}'.format(side)]
        _remove_multiplier(knee.spatialTransform.get_translation1())
        _remove_multiplier(knee.spatialTransform.get_translation2())
        _remove_multiplier(knee.spatialTransform.get_translation3())

        # update knee_l joint
        if self.verbose:
            print('updating knee {} joint...'.format(side))

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

        if self.verbose:
            print('knee {} splines:'.format(side))
            print(knee_spline_xk)
            print(knee_spline_yk)

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

        if self.verbose:
            print(
                'knee flexion: {:5.2f}'.format(
                    knee_joint.coordSets['knee_angle_{}'.format(side)].defaultValue
                    )
                )
        
        # update mesh
        if self.verbose:
            print('updating visual geometry...')
            
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
        if self.verbose:
            print('\nCUSTOMISING FOOT {}'.format(side.upper()))

        if (side!='l') and (side!='r'):
            raise ValueError('Invalid side')

        tibfib = self.LL.models['tibiafibula-'+side]
        femur = self.LL.models['femur-'+side]

        # scale foot bodies and joints
        # sf = scaler.calc_whole_body_scale_factors(
        #         self.LL, self._unit_scaling,
        #         )
        scaler.scale_body_mass_inertia(
            self.osimmodel.bodies['talus_{}'.format(side)],
            self._get_osimbody_scale_factors('talus_{}'.format(side))
            )
        scaler.scale_body_mass_inertia(
            self.osimmodel.bodies['calcn_{}'.format(side)],
            self._get_osimbody_scale_factors('calcn_{}'.format(side))
            )
        scaler.scale_body_mass_inertia(
            self.osimmodel.bodies['toes_{}'.format(side)],
            self._get_osimbody_scale_factors('toes_{}'.format(side))
            )
        
        scaler.scale_joint(
            self.osimmodel.joints['subtalar_{}'.format(side)],
            [
                self._get_osimbody_scale_factors('talus_{}'.format(side)),
                self._get_osimbody_scale_factors('calcn_{}'.format(side)),
            ],
            ['talus_{}'.format(side), 'calcn_{}'.format(side)]
            )
        scaler.scale_joint(
            self.osimmodel.joints['mtp_{}'.format(side)],
            [
                self._get_osimbody_scale_factors('calcn_{}'.format(side)),
                self._get_osimbody_scale_factors('toes_{}'.format(side)),
            ],
            ['calcn_{}'.format(side), 'toes_{}'.format(side)]
            )
        if self.verbose:
            print('scale factor: {}'.format(
                self._get_osimbody_scale_factors('talus_{}'.format(side))
                )
            )

        # remove multiplier functions from joint translations
        ankle = self.osimmodel.joints['ankle_{}'.format(side)]
        _remove_multiplier(ankle.spatialTransform.get_translation1())
        _remove_multiplier(ankle.spatialTransform.get_translation2())
        _remove_multiplier(ankle.spatialTransform.get_translation3())

        subtalar = self.osimmodel.joints['subtalar_{}'.format(side)]
        _remove_multiplier(subtalar.spatialTransform.get_translation1())
        _remove_multiplier(subtalar.spatialTransform.get_translation2())
        _remove_multiplier(subtalar.spatialTransform.get_translation3())

        mtp = self.osimmodel.joints['mtp_{}'.format(side)]
        _remove_multiplier(mtp.spatialTransform.get_translation1())
        _remove_multiplier(mtp.spatialTransform.get_translation2())
        _remove_multiplier(mtp.spatialTransform.get_translation3())

        # set ankle joint parent location in custom tibiafibula
        if self.verbose:
            print('updating ankle {} joint...'.format(side))

        ankle_centre = 0.5*(
            tibfib.landmarks['tibiafibula-MM'] + tibfib.landmarks['tibiafibula-LM']
            )
        self.osimmodel.joints['ankle_{}'.format(side)].locationInParent = \
            (tibfib.acs.map_local(ankle_centre[np.newaxis]).squeeze()*self._unit_scaling)+\
            self.ankle_offset

        if self.verbose:
            print(
                'location in parent: {}'.format(
                    self.osimmodel.joints['ankle_{}'.format(side)].locationInParent
                    )
                )

    def cust_osim_torso(self):
        if self.verbose:
            print('\nCUSTOMISING TORSO')

        pelvis = self.LL.models['pelvis']


        # scale torso inertial
        # sf = scaler.calc_whole_body_scale_factors(
        #         self.LL, self._unit_scaling,
        #         )
        sf = self._get_osimbody_scale_factors('torso')
        scaler.scale_body_mass_inertia(
            self.osimmodel.bodies['torso'], sf
            )
        if self.verbose:
            print('scale factor: {}'.format(sf))

        # remove multiplier functions from joint translations
        back = self.osimmodel.joints['back']
        _remove_multiplier(back.spatialTransform.get_translation1())
        _remove_multiplier(back.spatialTransform.get_translation2())
        _remove_multiplier(back.spatialTransform.get_translation3())

        # set back joint parent location in custom pelvis
        if self.verbose:
            print('updating back joint...')

        sacrum_top = pelvis.landmarks['pelvis-SacPlat']
        self.osimmodel.joints['back'].locationInParent = \
            (pelvis.acs.map_local(sacrum_top[np.newaxis]).squeeze()*self._unit_scaling)+\
            self.back_offset

        if self.verbose:
            print(
                'location in parent: {}'.format(
                    self.osimmodel.joints['back'].locationInParent
                    )
                )

    def write_cust_osim_model(self):
        self.osimmodel.save(
            os.path.join(str(self.config['osim_output_dir']), OSIM_FILENAME)
            )

    def customise(self):
        
        # model_scale_factors = self.scale_model()
        # self.scale_all_bodies()
        # self.recover_simmsplines()
        # model_scale_factors = self._calc_body_scale_factors()

        # for debugging: get original muscle optimal fibre lengths and tendon
        # slack lengths
        init_muscle_ofl = dict([(m.name, m.optimalFiberLength) for m in self.osimmodel.muscles.values()])
        init_muscle_tsl = dict([(m.name, m.tendonSlackLength) for m in self.osimmodel.muscles.values()])
        
        # prescale muscles to save their unscaled lengths
        self.prescale_muscles()

        # for debugging: get pre-scaled muscle optimal fibre lengths and tendon
        # slack lengths. Should not have changed
        prescale_muscle_ofl = dict([(m.name, m.optimalFiberLength) for m in self.osimmodel.muscles.values()])
        prescale_muscle_tsl = dict([(m.name, m.tendonSlackLength) for m in self.osimmodel.muscles.values()])

        # for m in self.osimmodel.muscles.values():
        #     print('{} tsl {} ofl {}'.format(m.name, m.tendonSlackLength, m.optimalFiberLength))

        # scale and modify bodies and joints
        self.cust_osim_pelvis()
        self.cust_osim_femur_l()
        self.cust_osim_femur_r()
        self.cust_osim_tibiafibula_l()
        self.cust_osim_tibiafibula_r()
        self.cust_osim_ankle_l()
        self.cust_osim_ankle_r()
        self.cust_osim_torso()

        # normalise the mass of each body against total subject mass (if provided)
        self.normalise_mass()

        # post-scale muscles to calculate their scaled lengths
        self.postscale_muscles()

        # for debugging: get scaled muscle optimal fibre lengths and tendon
        # slack lengths
        postscale_muscle_ofl = dict([(m.name, m.optimalFiberLength) for m in self.osimmodel.muscles.values()])
        postscale_muscle_tsl = dict([(m.name, m.tendonSlackLength) for m in self.osimmodel.muscles.values()])

        # for debugging: print out OFL and TSL changes through scaling
        if self.verbose:
            print('\nSCALED MUSCLE FIBRE PROPERTIES')
            for mn in sorted(self.osimmodel.muscles.keys()):
                print('{} OFL: {:8.6f} -> {:8.6f} -> {:8.6f}'.format(
                    mn, 
                    init_muscle_ofl[mn],
                    prescale_muscle_ofl[mn],
                    postscale_muscle_ofl[mn]
                    )
                )
            for mn in sorted(self.osimmodel.muscles.keys()):
                print('{} TSL: {:8.6f} -> {:8.6f} -> {:8.6f}'.format(
                    mn, 
                    init_muscle_tsl[mn],
                    prescale_muscle_tsl[mn],
                    postscale_muscle_tsl[mn]
                    )
                )

        # scale default markerset and add to model
        self.add_markerset()

        # write .osim file
        if self.config['write_osim_file']:
            self.write_cust_osim_model()

    def prescale_muscles(self):
        """
        Apply prescaling and scaling to muscles before bodies and joints are
        customised
        """
        state_0 = self.osimmodel._model.initSystem()
        scale_factors = scaler.calc_scale_factors_all_bodies(
            self.LL, self._unit_scaling, self.config['scale_other_bodies']
            )
        for m in self.osimmodel.muscles.values():
            m.preScale(state_0, *scale_factors)
            m.scale(state_0, *scale_factors)

    def postscale_muscles(self):
        """
        Postscale muscles after bodies and joints are customised to update
        optimal fiber lengths and tendon slack lengths
        """
        state_1 = self.osimmodel._model.initSystem()
        scale_factors = scaler.calc_scale_factors_all_bodies(
            self.LL, self._unit_scaling, self.config['scale_other_bodies']
            )
        for m in self.osimmodel.muscles.values():
            m.postScale(state_1, *scale_factors)

    def scale_model(self):
        model_sfs = scaler.calc_scale_factors_all_bodies(
            self.LL, self._unit_scaling, self.config['scale_other_bodies']
            )
        self.osimmodel.scale(self._osimmodel_init_state, *model_sfs)
        return model_sfs
    
    def add_markerset(self):
        """
        Add the default 2392 markerset to the customised osim model
        with customised marker positions.

        Markers in config['adj_marker_pairs'].keys() are placed in their
        corresponding input markers position.

        Else markers with bony landmark equivalents on the fieldwork model
        are assign the model landmarks with offset.

        Markers not matching the two above criteria are scaled according their
        body's scale factors.
        """
        vm = scaler.virtualmarker
        g2392_markers = vm._load_virtual_markers()[0]
        # maps of opensim names to fw names for markers and bodies
        osim2fw_markernames = dict([(it[1], it[0]) for it in vm.marker_name_map.items()])
        osim2fw_bodynames = dict([(it[1], it[0]) for it in OSIM_BODY_NAME_MAP.items()])

        adj_marker_pairs = self.config['adj_marker_pairs']
        if adj_marker_pairs is None:
            adj_marker_pairs = {}
        adj_model_marker_names = set(list(adj_marker_pairs.keys()))
        print('adj model markers:')
        for mm, mi in adj_marker_pairs.items():
            print('{} : {}'.format(mm,mi))

        def _local_coords(bodyname, landmarkname, global_coord=None, apply_offset=True):
            """
            Returns the local coordinates of a landmark
            """
            if global_coord is None:
                if landmarkname[-2:] in ('-l', '-r'):
                    _landmarkname = landmarkname[:-2]
                else:
                    _landmarkname = landmarkname
                global_coord = self.LL.models[bodyname].landmarks[_landmarkname]

            local_coords = self.LL.models[bodyname].acs.map_local(
                global_coord[np.newaxis,:]
                ).squeeze()

            if apply_offset:
                return self._unit_scaling*(local_coords + vm.marker_offsets[landmarkname])
            else:
                return self._unit_scaling*local_coords

        def _scale_marker(marker):
            """
            Scales the default opensim marker position by the scaling factor
            for its body
            """
            body_sf = self._get_osimbody_scale_factors(marker.bodyName)
            return marker.offset*body_sf

        self.markerset = opensim.MarkerSet()
        for osim_marker_name, marker0 in g2392_markers.items():
            new_offset = None

            # if define, adjust marker position to input marker
            if osim_marker_name in adj_model_marker_names:
                # move marker to input marker coordinates
                fw_body_name = osim2fw_bodynames[marker0.bodyName]
                input_marker_name = adj_marker_pairs[osim_marker_name]
                input_marker_coords = self.input_markers.get(input_marker_name)
                if input_marker_coords is None:
                    print(
                        'WARNING: {} not found in input markers. {} will not be adjusted.'.format(
                            input_marker_name, osim_marker_name
                            )
                        )
                else:
                    new_offset = _local_coords(
                        fw_body_name,
                        None,
                        global_coord=input_marker_coords,
                        apply_offset=False
                        )
            
            # if new marker position has not been defined by adjustment, then either set
            # as bony landmark coord or scale
            if new_offset is None:
                if osim_marker_name in osim2fw_markernames:
                    # if maker has fw equivalent move marker to fw landmark position with offset
                    fw_body_name = osim2fw_bodynames[marker0.bodyName]
                    fw_landmark_name = osim2fw_markernames[osim_marker_name]
                    new_offset = _local_coords(
                        fw_body_name,
                        fw_landmark_name,
                        apply_offset=True,
                        )
                else:
                    # else scale default
                    new_offset = _scale_marker(marker0)

            new_marker = osim.Marker(
                bodyname=marker0.bodyName,
                offset=new_offset
                )
            new_marker.name = marker0.name
            self.markerset.adoptAndAppend(new_marker._osimMarker)

        self.osimmodel._model.replaceMarkerSet(self._osimmodel_init_state, self.markerset)

    def normalise_mass(self):
        """
        Normalises the mass of each body so that total mass equals the given total mass.
        Doesn't do anything if config['subject_mass'] is None.
        """

        if self.config.get('subject_mass') is None:
            return

        if self.verbose:
            print('\nNORMALISING BODY MASSES')

        # if perserving reference model mass distribution, simply calculate a
        # uniform scaling factor for all bodies from original and target
        # subject mass
        if self.config['preserve_mass_distribution'] is True:
            if self.verbose:
                print('Preserving mass distribution')
            total_mass_0 = np.sum(self._original_segment_masses.values())
            target_mass = float(self.config['subject_mass'])
            mass_scaling = target_mass/total_mass_0
            for bname in self._original_segment_masses:
                b = self.osimmodel.bodies[bname]
                b.mass = self._original_segment_masses[bname]*mass_scaling
                if self.verbose:
                    print('{}: {:5.2f} kg'.format(b.name, b.mass))
        else:
            # calculate scaling factors for each body
            target_mass = float(self.config['subject_mass'])
            total_mass_0 = np.sum([float(b.mass) for b in self.osimmodel.bodies.values()])
            mass_scaling = target_mass/total_mass_0

            # scale mass for each body
            for b in self.osimmodel.bodies.values():
                b.mass = b.mass*mass_scaling
                if self.verbose:
                    print('{}: {:5.2f} kg'.format(b.name, b.mass))

        total_mass_1 = np.sum([float(b.mass) for b in self.osimmodel.bodies.values()])

        if self.verbose:
            print('Target Mass: {} kg'.format(target_mass))
            print('Unnormalised Mass: {} kg'.format(total_mass_0))
            print('Normalised Mass: {} kg'.format(total_mass_1))


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

def _remove_multiplier(owner):
    """
    Replace a components MultiplierFunction with the original function
    found in the MultiplierFunction instance.
    """
    newfunc = owner.getFunction()
    if newfunc.getConcreteClassName()=='MultiplierFunction':
        oldfunc = opensim.MultiplierFunction_safeDownCast(newfunc).getFunction()
        owner.setFunction(oldfunc.clone())

