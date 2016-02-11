"""
Tests for the gait2392geomcustomiser.py module
"""
import sys
import numpy as np
import copy
sys.path.append('../')
import gait2392geomcustomiser as g23
from gias2.fieldwork.field import geometric_field
from gias2.musculoskeletal.bonemodels import bonemodels

def _outputModelDict(LL):
    outputModelDict = dict([(m[0], m[1].gf) for m in LL.models.items()])
    return outputModelDict


# generate a custom left lower limb geometry
LL = bonemodels.LowerLimbLeftAtlas('lower_limb_left')
LL.load_bones()
LL.update_all_models([1.0,],[0,], [0,0,0,0,0,0], [0.0,0,0],[-0.0])
LLOutputModelDict = _outputModelDict(LL)

# test config file
config = {'osim_output_dir': './output/',  
          }

# instantiate customiser
cust = g23.Gait2392GeomCustomiser(config)
cust.set_model_gfields(LLOutputModelDict)

# customise each bone
cust.cust_osim_pelvis()
cust.cust_osim_femur_left()
cust.cust_osim_tibiafibula_left()


