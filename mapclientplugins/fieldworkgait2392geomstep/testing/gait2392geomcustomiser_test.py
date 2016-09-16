"""
Tests for the gait2392geomcustomiser.py module
"""
import os
import sys
sys.path.append('../')
import numpy as np
import copy
from gias2.fieldwork.field import geometric_field
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlas

import gait2392geomcustomiser as g23
reload(g23)

from lltransform import LLTransformData

SELF_DIRECTORY = os.path.split(__file__)[0]
_shapeModelFilenameRight = os.path.join(SELF_DIRECTORY, 'data/shape_models/LLP26_right_mirrored_from_left_rigid.pc')
_boneModelFilenamesRight = {'pelvis': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof'),
                                  os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.ens'),
                                  os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.mesh'),
                                  ),
                       'femur': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                 os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_quartic_flat.ens'),
                                 os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_quartic_flat.mesh'),
                                 ),
                       'patella': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                   os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_right.ens'),
                                   os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_right.mesh'),
                                   ),
                       'tibiafibula': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_cubic_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_right_cubic_flat.ens'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_right_cubic_flat.mesh'),
                                       ),
                        }
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


# generate a custom left lower limb geometry
ll_params = (
    [3.0,], # pc weights
    [0,], # pc modes
    [0,0,0,0,0,0], # pelvis rigid
    [0.5,0,0], # hip rot l
    [-0.5,0,0], # hip rot r
    [-0.2], # knee rot l
    [-1.0], # knee rot r
)

# LL = bonemodels.LowerLimbLeftAtlas('lower_limb_left')
LL = lowerlimbatlas.LowerLimbAtlas('lower_limb')
LL.ll_l.bone_files = _boneModelFilenamesLeft
LL.ll_l.combined_pcs_filename = _shapeModelFilenameLeft
LL.ll_r.bone_files = _boneModelFilenamesRight
LL.ll_r.combined_pcs_filename = _shapeModelFilenameRight
LL.load_bones()

LL.update_all_models(*ll_params)

# inputModelDict = _outputModelDict(LL)
# inputModelDict['femur-left'] = inputModelDict['femur']
# inputModelDict['tibiafibula-left'] = inputModelDict['tibiafibula']

# llt = LLTransformData()
# llt.pelvisRigid = ll_params[2]
# llt.hipRot = ll_params[3]
# llt.kneeRot = ll_params[4]

# test config file
output_dir = str(os.path.join(os.path.split(__file__)[0], 'output/'))
config = {'osim_output_dir': output_dir,
          'in_unit': 'mm',
          'out_unit': 'm',
          'scale_other_bodies': True,
          'write_osim_file': True,
          'subject_mass': 80.0,
          'preserve_mass_distribution': True,
          'adj_marker_pairs': None,
          }

# instantiate customiser
cust = g23.Gait2392GeomCustomiser(config, ll=LL)
cust.init_osim_model()
# cust.set_lowerlimb_gfields(inputModelDict)
# cust.ll_transform = llt

# customise each bone
# cust.cust_osim_pelvis()
# cust.cust_osim_femur_left()
# cust.cust_osim_tibiafibula_left()
# cust.cust_osim_ankle_left()
cust.customise()

# print('writing')
# write out customised osim file
cust.write_cust_osim_model()


# knee_angles = g23.calc_knee_angles(cust.LL)
# print(knee_angles)

# view customised opensim model
# omodel = cust.osimmodel
# v = omodel.view_init_state()