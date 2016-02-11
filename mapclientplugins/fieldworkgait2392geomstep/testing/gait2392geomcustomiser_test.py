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
reload(g23)

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


# generate a custom left lower limb geometry
LL = bonemodels.LowerLimbLeftAtlas('lower_limb_left')
LL.bone_files = _boneModelFilenamesLeft
LL.combined_pcs_filename = _shapeModelFilenameLeft
LL.load_bones()
LL.update_all_models([1.0,],[0,], [0,0,0,0,0,0], [0.0,0,0],[-0.0])
LLOutputModelDict = _outputModelDict(LL)

# test config file
output_dir = str(os.path.join(os.path.split(__file__)[0], 'output/'))
config = {'osim_output_dir': output_dir,  
          }

# instantiate customiser
cust = g23.Gait2392GeomCustomiser(config)
cust.set_model_gfields(LLOutputModelDict)

# customise each bone
cust.cust_osim_pelvis()
cust.cust_osim_femur_left()
cust.cust_osim_tibiafibula_left()


