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

OpenSim Gait2392 bodies customisation
"""

import os
import opensim
import numpy as np
import copy

from gias2.mesh import vtktools
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlas
from gias2.musculoskeletal import osim

# from mapclientplugins.fieldworkgait2392geomstep import scaler
from . import scaler

# =============================================================================#
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
                        'LH': range(73, 146),
                        'sac': range(146, 260),
                        }
PELVIS_BASISTYPES = {'tri10': 'simplex_L3_L3', 'quad44': 'quad_L3_L3'}
TIBFIB_SUBMESHES = ('tibia', 'fibula')
TIBFIB_SUBMESH_ELEMS = {'tibia': range(0, 46),
                        'fibula': range(46, 88),
                        }
TIBFIB_BASISTYPES = {'tri10': 'simplex_L3_L3', 'quad44': 'quad_L3_L3'}

GEOM_DIR = '../Geometry'

VALID_UNITS = ('nm', 'um', 'mm', 'cm', 'm', 'km')
VALID_MODEL_MARKERS = sorted(list(scaler.virtualmarker.markers.keys()))


# =============================================================================#
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
        'm': 1.0,
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

    return unit_vals[in_unit] / unit_vals[out_unit]


# Opensim coordinate systems for bodies. We only seem to be doing this for the
# femur and tibia-fibula bones. Do we need to for the other body components?
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


def _calc_knee_spline_coords(ll, flex_angles):
    """
    Calculates the cubic spline values for the knee joint through specified
    angles. The values are the coordinates of the tibia frame origin relative
    to the femur frame.

    inputs
    ======
    ll : LowerLimbLeftAtlas instance
    flex_angles : 1d ndarray
        A list of n knee angles at which to sample tibia location relative to
        the femur ACS. Only flexion supported in 2392.

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
        _ll.update_tibiafibula([a, 0, 0])
        tib_o = 0.5 * (_ll.models['tibiafibula'].landmarks['tibiafibula-LC'] +
                       _ll.models['tibiafibula'].landmarks['tibiafibula-MC']
                       )
        tib_os.append(tib_o)

    update_femur_opensim_acs(_ll.models['femur'])
    y = _ll.models['femur'].acs.map_local(np.array(tib_os))
    # y = np.array([y[:,2], y[:,1], y[:,0]]).T # reverse dims
    return y


# ========================================================================== #
class Gait2392GeomCustomiser(object):
    gfield_disc = (6, 6)
    ankle_offset = np.array([0., -0.01, 0.])
    back_offset = np.array([0., 0.01, 0.])

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
        self.uniform_scaling = 1.0
        self.pelvis_scaling = 1.0
        self.femur_scaling = 1.0
        self.petalla_scaling = 1.0
        self.tibfib_scaling = 1.0
        self.LL = None  # Gias2 LowerLimbAtlas object.
        self._hasInputLL = False
        self._workflow_location = None
        self.osimmodel = None  # OpenSim Model.
        self.markerset = None  # MarkerSet associated with OpenSim Model.
        self.input_markers = {}  # Input Marker name : input Marker coords.
        self.verbose = verbose
        self._unit_scaling = dim_unit_scaling(
            self.config['in_unit'], self.config['out_unit']
        )

        if gfieldsdict is not None:
            self.set_lowerlimb_gfields(gfieldsdict)
        if ll is not None:
            self.set_lowerlimb_atlas(ll)

        # These hold the relevant scaling information that will be used to
        # scale the Model. They are created when the LowerLimbAtlas is set.
        self.scale_set = []
        self.scale_factors = {}

        self._osimmodel_state = None
        self._original_segment_masses = None

    def init_osim_model(self):
        self.osimmodel = osim.Model(TEMPLATE_OSIM_PATH)

        self._osimmodel_state = self.osimmodel.init_system()

        self._original_segment_masses = dict([
            (b.name, b.mass) for b in self.osimmodel.bodies.values()])

    def set_workflow_location(self, location):
        self._workflow_location = location

    def _check_geom_path(self):
        """
        Check that the directory for geom meshes exists. If not, create it.
        """
        output_directory = self.config['osim_output_dir']
        if not os.path.isabs(output_directory):
            output_directory = os.path.join(self._workflow_location, output_directory)

        geom_dir = os.path.join(output_directory, GEOM_DIR)

        if not os.path.isdir(geom_dir):
            os.mkdir(geom_dir)

    def set_lowerlimb_atlas(self, ll):
        """
        Set the Gias2 LowerLimbAtlas and calculates the set of Scales (and from
        this the set of scale factors) that will be used to scale the model.
        """
        self.LL = ll
        self._hasInputLL = True

        update_femur_opensim_acs(self.LL.ll_l.models['femur'])
        update_tibiafibula_opensim_acs(self.LL.ll_l.models['tibiafibula'])
        update_femur_opensim_acs(self.LL.ll_r.models['femur'])
        update_tibiafibula_opensim_acs(self.LL.ll_r.models['tibiafibula'])

        self.scale_set = scaler.calc_scale_factors_all_bodies(
            ll, self._unit_scaling)
        for scale in self.scale_set:
            scale_factors = scale.scaleFactors
            self.scale_factors[scale.segmentName] = scale_factors

    def set_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """
        self.set_2side_lowerlimb_gfields(gfieldsdict)

    # This method looks like it could be condensed significantly.
    def set_2side_lowerlimb_gfields(self, gfieldsdict):
        """
        Instantiate the lower limb object using input models
        """

        print('Setting lower limb fieldwork models')
        print(gfieldsdict.keys())

        # left
        if not self._hasInputLL:
            print('creating new left lowerlimb model')
            ll_l = bonemodels.LowerLimbLeftAtlas('left lower limb')
            ll_l.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            ll_l.set_bone_gfield('femur', gfieldsdict['femur-l'])
            ll_l.set_bone_gfield('patella', gfieldsdict['patella-l'])
            ll_l.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-l'])
        else:
            ll_l = self.LL.ll_l
            if 'pelvis' in gfieldsdict:
                print('setting custom pelvis (l)')
                ll_l.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            if 'femur-l' in gfieldsdict:
                print('setting custom femur-l')
                ll_l.set_bone_gfield('femur', gfieldsdict['femur-l'])
            if 'patella-l' in gfieldsdict:
                print('setting custom patella-l')
                ll_l.set_bone_gfield('patella', gfieldsdict['patella-l'])
            if 'tibiafibula-l' in gfieldsdict:
                print('setting custom tibiafibula-l')
                ll_l.set_bone_gfield(
                    'tibiafibula', gfieldsdict['tibiafibula-l'])

        update_femur_opensim_acs(ll_l.models['femur'])
        update_tibiafibula_opensim_acs(ll_l.models['tibiafibula'])

        # right
        if not self._hasInputLL:
            print('creating new right lowerlimb model')
            ll_r = bonemodels.LowerLimbLeftAtlas('right lower limb')
            ll_r.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            ll_r.set_bone_gfield('femur', gfieldsdict['femur-r'])
            ll_r.set_bone_gfield('patella', gfieldsdict['patella-r'])
            ll_r.set_bone_gfield('tibiafibula', gfieldsdict['tibiafibula-r'])
        else:
            ll_r = self.LL.ll_r
            if 'pelvis' in gfieldsdict:
                print('setting custom pelvis (r)')
                ll_r.set_bone_gfield('pelvis', gfieldsdict['pelvis'])
            if 'femur-r' in gfieldsdict:
                print('setting custom femur-r')
                ll_r.set_bone_gfield('femur', gfieldsdict['femur-r'])
            if 'patella-r' in gfieldsdict:
                print('setting custom patella-r')
                ll_r.set_bone_gfield('patella', gfieldsdict['patella-r'])
            if 'tibiafibula-r' in gfieldsdict:
                print('setting custom tibiafibula-r')
                ll_r.set_bone_gfield(
                    'tibiafibula', gfieldsdict['tibiafibula-r'])

        update_femur_opensim_acs(ll_r.models['femur'])
        update_tibiafibula_opensim_acs(ll_r.models['tibiafibula'])

        # 2side
        if not self._hasInputLL:
            self.LL = lowerlimbatlas.LowerLimbAtlas('lower limb')
            self.LL.ll_l = ll_l
            self.LL.ll_r = ll_r

        self.LL.update_model_dict()

    def _save_stl(self, gf, filename, bodycoordmapper):
        v, f = gf.triangulate(self.gfield_disc)
        v_local = bodycoordmapper(v)
        v_local *= self._unit_scaling
        vtkwriter = vtktools.Writer(
            v=v_local,
            f=f,
            filename=filename,
        )
        vtkwriter.writeSTL()

    def _get_osim_knee_spline_xk(self, side):
        """
        Get the SimmSpline x values from the translation functions of the
        gati2392 knee
        """
        if (side != 'l') and (side != 'r'):
            raise ValueError('Invalid side')

        if side == 'l':
            kj = self.osimmodel.joints['knee_l']
        else:
            kj = self.osimmodel.joints['knee_r']

        t1x = kj.getSimmSplineParams('translation1')[0]
        t2x = kj.getSimmSplineParams('translation2')[0]
        return t1x, t2x

    def _set_osim_knee_spline_xyk(self, x, y, side):
        if (side != 'l') and (side != 'r'):
            raise ValueError('Invalid side')

        if side == 'l':
            kj = self.osimmodel.joints['knee_l']
        else:
            kj = self.osimmodel.joints['knee_r']

        kj.updateSimmSplineParams('translation1', x[0], y[0])
        kj.updateSimmSplineParams('translation2', x[1], y[1])

    def write_cust_osim_model(self):
        """
        Save the current OpenSim::Model object to an .osim XML file. The output
        directory is given by the plugin configuration.
        """
        output_directory = self.config['osim_output_dir']
        if not os.path.isabs(output_directory):
            output_directory = os.path.join(self._workflow_location, output_directory)

        self.osimmodel.save(os.path.join(output_directory, OSIM_FILENAME))

    # We may need to remove the scale factors that have been applied to the
    # Model before (in) this method. But only for the "generated" mesh files.
    def save_mesh_files(self):
        """
        Save the updated mesh files associated with the Model.
        """

        self._check_geom_path()

        name_list = ["femur-l", "tibiafibula-l", "femur-r", "tibiafibula-r",
                     "pelvis"]
        for body_name in name_list:
            mesh_list = self.save_mesh_file(body_name)
            osim_body = self.osimmodel.bodies[OSIM_BODY_NAME_MAP[body_name]]

            osim_body.setDisplayGeometryFileName(mesh_list)

    def save_mesh_file(self, body_name):

        def _split_tibia_fibula_gfs(tib_fib_gf):
            tib = tib_fib_gf.makeGFFromElements(
                'tibia',
                TIBFIB_SUBMESH_ELEMS['tibia'],
                TIBFIB_BASISTYPES,
            )
            fib = tib_fib_gf.makeGFFromElements(
                'fibula',
                TIBFIB_SUBMESH_ELEMS['fibula'],
                TIBFIB_BASISTYPES,
            )

            return tib, fib

        def _split_pelvis_gfs(pelvis_gf):
            lhgf = pelvis_gf.makeGFFromElements(
                'hemipelvis-left',
                PELVIS_SUBMESH_ELEMS['LH'],
                PELVIS_BASISTYPES
            )
            sacgf = pelvis_gf.makeGFFromElements(
                'sacrum',
                PELVIS_SUBMESH_ELEMS['sac'],
                PELVIS_BASISTYPES
            )
            rhgf = pelvis_gf.makeGFFromElements(
                'hemipelvis-right',
                PELVIS_SUBMESH_ELEMS['RH'],
                PELVIS_BASISTYPES
            )
            return lhgf, sacgf, rhgf

        if body_name == "pelvis":
            file_names = ["sacrum.stl", "r_pelvis.stl", "l_pelvis.stl"]
            gf_list = _split_pelvis_gfs(self.LL.models['pelvis'].gf)
            acs_map_local = self.LL.models['pelvis'].acs.map_local

        elif body_name[0:5] == "femur":
            file_names = [body_name[-1] + "_femur.stl"]
            ll_model = self.LL.models[body_name]
            gf_list = [ll_model.gf]
            acs_map_local = ll_model.acs.map_local

        else:
            file_names = [body_name[-1] + "_tibia.stl",
                          body_name[-1] + "_fibula.stl"]
            ll_model = self.LL.models[body_name]
            gf_list = _split_tibia_fibula_gfs(ll_model.gf)
            acs_map_local = ll_model.acs.map_local

        output_directory = self.config['osim_output_dir']
        if not os.path.isabs(output_directory):
            output_directory = os.path.join(self._workflow_location, output_directory)

        for i in range(len(file_names)):
            stl_full_path = os.path.join(output_directory, GEOM_DIR, file_names[i])

            self._save_stl(gf_list[i], stl_full_path, acs_map_local)

        return file_names

    def update_knee_splines(self):
        self._update_side_knee_splines('l')
        self._update_side_knee_splines('r')

    def _update_side_knee_splines(self, side):
        """
        Update knee joint splines to reflect the changes to the Model.
        """

        if side == 'l':
            knee_spline_xk_1, knee_spline_xk_2 = self._get_osim_knee_spline_xk(
                'l')
            ll = self.LL.ll_l
        else:
            knee_spline_xk_1, knee_spline_xk_2 = self._get_osim_knee_spline_xk(
                'r')
            ll = self.LL.ll_r

        knee_spline_xk = [knee_spline_xk_1, knee_spline_xk_2]

        # Evaluate tib coord at xks.
        knee_spline_yk_1 = _calc_knee_spline_coords(
            ll, knee_spline_xk_1) * self._unit_scaling
        knee_spline_yk_2 = _calc_knee_spline_coords(
            ll, knee_spline_xk_2) * self._unit_scaling

        knee_spline_yk = [knee_spline_yk_1[:, 0], knee_spline_yk_2[:, 1]]
        # Set new spline yks.
        self._set_osim_knee_spline_xyk(knee_spline_xk, knee_spline_yk, side)

        if self.verbose:
            print('knee {} splines:'.format(side))
            print(knee_spline_xk)
            print(knee_spline_yk)

    # Note: this doesn't remove any MultiplierFunction wrappers from the
    # Muscle MovingPathPoints.
    # Perhaps this method could be integrated into update_joints to save on
    # function calls.
    # Note: this function will remove any scale_factors, if these are
    # important they should be applied to the values explicitly (if possible)
    # before using this method.
    def remove_multipliers(self):
        """
        Replace each components MultiplierFunction with the original function
        found in the MultiplierFunction instance. These wrappers seem to have
        been added to translational axis SimmSpline functions during scaling.
        """

        def _remove_multiplier(owner):
            current_func = owner.get_function()
            if current_func.getConcreteClassName() == 'MultiplierFunction':
                m_func = opensim.MultiplierFunction.safeDownCast(
                    current_func)
                new_func = m_func.getFunction()
                owner.set_function(new_func.clone())

        for joint_name, joint in self.osimmodel.joints.items():
            if joint_name == 'ground_pelvis':
                continue
            _remove_multiplier(joint.spatialTransform.get_translation1())
            _remove_multiplier(joint.spatialTransform.get_translation2())
            _remove_multiplier(joint.spatialTransform.get_translation3())

    def customise(self):
        """
        Customise the Model. See the plugin README for more information.
        """

        # Scale the Model.
        self.scale_model()

        # Update each joint according to its position in the LowerLimb model.
        self.remove_multipliers()
        self.update_joints()
        self.update_knee_splines()

        # Save the geometry mesh files for the updated Model.
        self.save_mesh_files()

        # Scale the default MarkerSet and add it to the OpenSim Model.
        self.add_markerset()

        # Save the Model to an .osim XML file.
        if self.config['write_osim_file']:
            self.write_cust_osim_model()

    def scale_model(self):
        """
        This method handles most of the scaling for the model. It returns a
        list of the Scale objects that were used to scale the Model.
        """
        # We need to confirm that the body masses are correctly scaled for
        # the gait2392 model when the 'subject_mass' is supplied. If not, look
        # at the normalise_mass method (removed).
        if self.config['subject_mass'] is None:
            mass = -1.0
        else:
            mass = self.config['subject_mass']

        self.osimmodel.scale(self._osimmodel_state, *self.scale_set,
                             preserve_mass_distribution=self.config[
                                 'preserve_mass_distribution'], mass=mass)

        self._osimmodel_state = self.osimmodel.init_system()

        return self.scale_set

    def update_joints(self):
        """
        Location and LocationInParent are modified according to the position of
        the joints in the Gias2 model.
        """
        # Only the following Joints have Gias2 associates. All other joints
        # should be scaled by the relevant body scale factors (scale_model).
        joint_names = ["ground_pelvis", "hip_l", "hip_r", "knee_l", "knee_r",
                       "ankle_l", "ankle_r", "back"]
        pelvis = self.LL.models['pelvis']
        location = None
        offset = np.array([0., 0., 0.])
        for joint_name in joint_names:
            joint = self.osimmodel.joints[joint_name]

            if joint_name == "ground_pelvis":
                location_in_parent = self.LL.models['pelvis'].acs.o
                location = np.array((0, 0, 0), dtype=float)

                tilt, _list, rot = self.LL.pelvis_rigid[3:]
                joint.coordSets['pelvis_tilt'].defaultValue = tilt
                joint.coordSets['pelvis_list'].defaultValue = _list
                joint.coordSets['pelvis_rotation'].defaultValue = rot

            elif joint_name == "hip_l" or joint_name == "hip_r":
                side = joint_name[-1]
                hjc = pelvis.landmarks['pelvis-' + side.upper() + 'HJC']
                femur = self.LL.models['femur-' + side]
                location_in_parent = pelvis.acs.map_local(hjc[np.newaxis])[0]
                location = femur.acs.map_local(hjc[np.newaxis])[0]

                if side == 'l':
                    flex, rot, add = self.LL.hip_rot_l
                else:
                    flex, rot, add = self.LL.hip_rot_r
                joint.coordSets[f'hip_flexion_{side}'].defaultValue = flex
                joint.coordSets[f'hip_adduction_{side}'].defaultValue = add
                joint.coordSets[f'hip_rotation_{side}'].defaultValue = rot

            elif joint_name == "knee_l" or joint_name == "knee_r":
                side = joint_name[-1]
                tibfib = self.LL.models['tibiafibula-' + side]
                femur = self.LL.models['femur-' + side]
                # In the outdated version of this method, the following three
                # values were being calculated but not used. Are they useful?
                kjc = 0.5 * (femur.landmarks['femur-MEC'] +
                             femur.landmarks['femur-LEC'])
                tpc = 0.5 * (tibfib.landmarks['tibiafibula-MC'] +
                             tibfib.landmarks['tibiafibula-LC'])
                _d = -np.sqrt(((kjc - tpc) ** 2.0).sum())
                location_in_parent = np.array([0, 0, 0], dtype=float)
                location = np.array([0, 0, 0], dtype=float)

                if side == 'l':
                    flex, add = self.LL.knee_rot_l
                else:
                    flex, add = self.LL.knee_rot_r
                joint.coordSets[f'knee_angle_{side}'].defaultValue = flex

            elif joint_name == "ankle_l" or joint_name == "ankle_r":
                side = joint_name[-1]
                tibfib = self.LL.models['tibiafibula-' + side]
                ankle_centre = 0.5 * (tibfib.landmarks['tibiafibula-MM'] +
                                      tibfib.landmarks['tibiafibula-LM'])
                location_in_parent = tibfib.acs.map_local(
                    ankle_centre[np.newaxis]).squeeze()
                offset = self.ankle_offset

                # Shouldn't we update the ankle angle Coordinate default value?
                # The LowerLimbAtlas class doesn't appear to have a method to
                # get this value from the model. Same for the 'back' Joint.

            else:
                sacrum_top = pelvis.landmarks['pelvis-SacPlat']
                location_in_parent = pelvis.acs.map_local(
                    sacrum_top[np.newaxis]).squeeze()
                offset = self.back_offset

            joint.locationInParent = (
                    location_in_parent * self._unit_scaling) + offset
            if location is not None:
                joint.location = location * self._unit_scaling

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
        g2392_markers = vm.g2392_markers

        # Maps of opensim names to fw names for markers and bodies.
        osim2fw_markernames = dict(
            [(it[1], it[0]) for it in vm.marker_name_map.items()])
        osim2fw_bodynames = dict(
            [(it[1], it[0]) for it in OSIM_BODY_NAME_MAP.items()])

        adj_marker_pairs = self.config['adj_marker_pairs']
        if adj_marker_pairs is None:
            adj_marker_pairs = {}
        adj_model_marker_names = set(list(adj_marker_pairs.keys()))
        print('adj model markers:')
        for mm, mi in adj_marker_pairs.items():
            print('{} : {}'.format(mm, mi))

        def _local_coords(bodyname, landmarkname, global_coord=None,
                          apply_offset=True):
            """
            Returns the local coordinates of a landmark
            """
            if global_coord is None:
                if landmarkname[-2:] in ('-l', '-r'):
                    _landmarkname = landmarkname[:-2]
                else:
                    _landmarkname = landmarkname
                global_coord = self.LL.models[bodyname].landmarks[
                    _landmarkname]

            local_coords = self.LL.models[bodyname].acs.map_local(
                global_coord[np.newaxis, :]
            ).squeeze()

            if apply_offset:
                return self._unit_scaling * (local_coords + vm.marker_offsets[
                    landmarkname])
            else:
                return self._unit_scaling * local_coords

        def _scale_marker(marker):
            """
            Scales the default opensim marker position by the scaling factor
            for its body
            """

            body_sf = self.scale_factors[marker.frame_name]
            return marker.location * body_sf

        self.markerset = opensim.MarkerSet()
        for osim_marker_name, marker0 in g2392_markers.items():
            new_offset = None

            # If defined, adjust marker position to input marker.
            if osim_marker_name in adj_model_marker_names:
                # Move marker to input marker coordinates.
                fw_body_name = osim2fw_bodynames[marker0.frame_name]
                input_marker_name = adj_marker_pairs[osim_marker_name]
                input_marker_coords = self.input_markers.get(input_marker_name)
                if input_marker_coords is None:
                    print(
                        'WARNING: {} not found in input markers. {} will not '
                        'be adjusted.'.format(
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

            # If new marker position has not been defined by adjustment, then
            # either set as bony landmark coord or scale.
            if new_offset is None:
                if osim_marker_name in osim2fw_markernames:
                    # If maker has fw equivalent move marker to fw landmark
                    # position with offset.
                    fw_body_name = osim2fw_bodynames[marker0.frame_name]
                    fw_landmark_name = osim2fw_markernames[osim_marker_name]
                    new_offset = _local_coords(
                        fw_body_name,
                        fw_landmark_name,
                        apply_offset=True,
                    )
                else:
                    # Else scale default.
                    new_offset = _scale_marker(marker0)

            # Pre-pend "/bodyset/" back on to Marker's frame_name for writing.
            full_frame_name = "/bodyset/" + marker0.frame_name

            new_marker = osim.Marker(
                name=marker0.name,
                frame_name=full_frame_name,
                location=new_offset
            )
            self.markerset.adoptAndAppend(new_marker.get_osim_marker())

        self.osimmodel.set_marker_set(self.markerset)
