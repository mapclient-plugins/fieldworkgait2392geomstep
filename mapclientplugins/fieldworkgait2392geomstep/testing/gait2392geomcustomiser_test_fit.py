"""
Tests for the gait2392geomcustomiser.py module
"""
import os
import sys
import numpy as np
import copy
sys.path.append('../')
import gait2392geomcustomiser as g23
from gias2.fieldwork.field import geometric_field
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfit
from gias2.musculoskeletal import mocap_landmark_preprocess
from gias2.visualisation import fieldvi
reload(g23)

from lltransform import LLTransformData
import trcdata

SELF_DIRECTORY = os.path.split(__file__)[0]
_shapeModelFilenameLeft = os.path.join(SELF_DIRECTORY, 'data/shape_models/LLP26_rigid.pc')
_boneModelFilenamesLeft = {'pelvis': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof'),
                                  os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.ens'),
                                  os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.mesh'),
                                  ),
                       'femur': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_mean_rigid_LLP26.geof'),
                                 os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_quartic_flat.ens'),
                                 os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_quartic_flat.mesh'),
                                 ),
                       'patella': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_left_mean_rigid_LLP26.geof'),
                                   os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_left.ens'),
                                   os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_left.mesh'),
                                   ),
                       'tibiafibula': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_cubic_left_mean_rigid_LLP26.geof'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_left_cubic_flat.ens'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_left_cubic_flat.mesh'),
                                       ),
                       }

def _outputModelDict(LL):
    outputModelDict = dict([(m[0], m[1].gf) for m in LL.models.items()])
    return outputModelDict

#===================================================#
# fitting Parameters                                #
#===================================================#
model_name = 'lower_limb_left_LLP26'

# mocap landmarks
load_mocap = True
mocap_file = 'data/example_static1.trc'
mocap_frame = 5
target_landmark_names = ['L.ASIS', 'R.ASIS', 'V.Sacral', 'L.Knee.Medial',
                        'L.Knee', 'L.Ankle.Medial', 'L.Ankle',
                        ]
source_landmark_names = ['pelvis-LASIS', 'pelvis-RASIS', 'pelvis-Sacral',
                         'femur-MEC', 'femur-LEC', 'tibiafibula-MM',
                         'tibiafibula-LM',
                        ]

# fitting parameters
pc_modes = [0,]
mweight = 0.1
min_args = {'method':'L-BFGS-B',
            'jac':False,
            'bounds':None, 'tol':1e-6,
            'options':{'eps':1e-5},
            }

#===================================================#
# initialise                                        #
#===================================================#
# load lower limb model
LL = bonemodels.LowerLimbLeftAtlas(model_name)
LL.bone_files = _boneModelFilenamesLeft
LL.combined_pcs_filename = _shapeModelFilenameLeft
LL.load_bones()

# load target landmarks
trc = trcdata.TRCData()
trc.load(mocap_file)
all_landmarks = trc.get_frame(mocap_frame)
target_landmark_coords_raw = np.array([all_landmarks[l] for l in target_landmark_names])
target_landmark_coords = np.array(mocap_landmark_preprocess.preprocess_lower_limb(
                                  5.0, 10.0, *target_landmark_coords_raw))


#===================================================#
# Fit                                               #
#===================================================#
# single stage fit
fitting_xs,\
opt_landmark_dist,\
opt_landmark_rmse,\
min_info = lowerlimbatlasfit.fit(
                LL, target_landmark_coords, source_landmark_names,
                pc_modes, mweight, minimise_args=min_args)
fitted_landmark_coords = min_info['opt_source_landmarks']
print('Fitted landmark RMSE: {}'.format(opt_landmark_rmse))


inputModelDict = _outputModelDict(LL)

# llt = LLTransformData()
# llt.pelvisRigid = ll_params[2]
# llt.hipRot = ll_params[3]
# llt.kneeRot = ll_params[4]

# test config file
output_dir = str(os.path.join(os.path.split(__file__)[0], 'output/'))
config = {'osim_output_dir': output_dir,
          'convert_mm_to_m': True,
          'scale_other_bodies': False,
          'write_osim_file': True,
          }

# instantiate customiser
cust = g23.Gait2392GeomCustomiser(config)
cust.set_left_lowerlimb_gfields(inputModelDict)
# cust.ll_transform = llt

# customise each bone
# cust.cust_osim_pelvis()
# cust.cust_osim_femur_left()
# cust.cust_osim_tibiafibula_left()
# cust.cust_osim_ankle_left()
cust.customise()

# print('writing')
# write out customised osim file
# cust.write_cust_osim_model()

knee_angles = g23.calc_knee_angles(cust.LL)
print(knee_angles)

#============================================#
# Visualisation                              #
#============================================#
# V = fieldvi.Fieldvi()
# V.GFD = [6,6]
# V.displayGFNodes = False
# for mn, m in LL.models.items():
#     mgfeval = geometric_field.makeGeometricFieldEvaluatorSparse(m.gf, V.GFD)
#     V.addGeometricField(mn, m.gf, mgfeval, V.GFD)

# V.addData('pelvis_acs', LL.models['pelvis'].acs.unit_array, scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('femur_acs', LL.models['femur'].acs.unit_array, scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('patella_acs', LL.models['patella'].acs.unit_array, scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('tibiafibula_acs', LL.models['tibiafibula'].acs.unit_array, scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('hip_acs', LL._get_hip_cs(), scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('knee_acs', LL._get_knee_cs(), scalar=(0,1,2,3), renderArgs={'scale_mode':'none'})
# V.addData('fitted landmarks', fitted_landmark_coords,
#           scalar=range(len(fitted_landmark_coords)),
#           renderArgs={'scale_mode':'none', 'scale_factor':6.0})
# V.addData('target landmarks', target_landmark_coords,
#           scalar=range(len(target_landmark_names)),
#           renderArgs={'scale_mode':'none', 'scale_factor':6.0})
# V.configure_traits()
# V.scene.background = (0,0,0)