"""
Module for reading Gait2392's virtual markerset
"""

import os
import opensim
import numpy as np

from gias2.musculoskeletal.osim import Marker

SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
MARKERSET_PATH = str(os.path.join(SELF_DIR, 'data', 'gait2392_Scale_MarkerSet.xml'))
MARKER_OFFSET_PATH = str(os.path.join(SELF_DIR, 'data/', 'marker_offsets.dat'))
try:
    opensim_version = getattr(opensim, '__version__')
except AttributeError:
    opensim_version = None

# dictionary mapping fieldwork landmark names to gait2392's virtual marker
# names
marker_name_map = {
    'pelvis-RASIS': 'R.ASIS',
    'pelvis-LASIS': 'L.ASIS',
    'pelvis-Sacral': 'V.Sacral',
    'femur-LEC-r': 'R.Knee.Lat',
    'femur-MEC-r': 'R.Knee.Med',
    'tibiafibula-LM-r': 'R.Ankle.Lat',
    'tibiafibula-MM-r': 'R.Ankle.Med',
    'femur-LEC-l': 'L.Knee.Lat',
    'femur-MEC-l': 'L.Knee.Med',
    'tibiafibula-LM-l': 'L.Ankle.Lat',
    'tibiafibula-MM-l': 'L.Ankle.Med',
    'pelvis-LHJC': 'L.HJC',
    'pelvis-RHJC': 'R.HJC',
    'femur-HC-l': 'L.FHC',
    'femur-HC-r': 'R.FHC',
    'tibiafibula-KJC-l': 'L.Tib.KJC',
    'tibiafibula-KJC-r': 'R.Tib.KJC',
}

def _load_virtual_markers():
    markers = {}
    marker_coords = {}

    if opensim_version==4.0:
        _dummy_model = opensim.Model()
        _osim_markerset = opensim.MarkerSet(_dummy_model, MARKERSET_PATH)
    else:
        _osim_markerset = opensim.MarkerSet(MARKERSET_PATH)

    for mi in range(_osim_markerset.getSize()):
        osim_marker = _osim_markerset.get(mi)

        # create a copy because the markerset and its marers only exists
        # in the function
        _v = opensim.Vec3()
        osim_marker.getOffset(_v)
        offset = np.array([_v.get(i) for i in range(3)])
        marker = Marker(
                    name=osim_marker.getName(),
                    bodyname=osim_marker.getBodyName(),
                    offset=offset,
                    )
        markers[marker.name] = marker
        marker_coords[marker.name] = marker.offset

    return markers, marker_coords

def _add_synthetic_markers(marker, marker_coords):
    # add hip joint centres to markers (joint coords taken from gait2392)
    markers['L.HJC'] = Marker(
        bodyname='pelvis',
        offset=(-0.0707, -0.0661, -0.0835)
        )
    marker_coords['L.HJC'] = markers['L.HJC'].offset
    markers['R.HJC'] = Marker(
        bodyname='pelvis',
        offset=(-0.0707, -0.0661, 0.0835)
        )
    marker_coords['R.HJC'] = markers['R.HJC'].offset

    # add femur head centre to markers
    markers['L.FHC'] = Marker(bodyname='femur_l', offset=(0,0,0))
    marker_coords['L.FHC'] = markers['L.FHC'].offset
    markers['R.FHC'] = Marker(bodyname='femur_r', offset=(0,0,0))
    marker_coords['R.FHC'] = markers['R.FHC'].offset

    # add knee centre in tibia frame to markers
    markers['L.Tib.KJC'] = Marker(bodyname='tibia_l', offset=(0,0,0))
    marker_coords['L.Tib.KJC'] = markers['L.Tib.KJC'].offset
    markers['R.Tib.KJC'] = Marker(bodyname='tibia_r', offset=(0,0,0))
    marker_coords['R.Tib.KJC'] = markers['R.Tib.KJC'].offset

def _load_marker_offsets():
    """
    Read from file the offsets between the opensim virtual markers and 
    their corresponding fieldwork landmarks
    """
    marker_offsets = {}
    with open(MARKER_OFFSET_PATH, 'r') as f:
        lines = f.readlines()

    for l in lines:
        if l[0]=='#':
            pass
        else:
            words = l.split()
            if len(words)==4:
                marker_offsets[words[0]] = np.array([float(x) for x in words[1:]])

    return marker_offsets

# load up virtual markers from file
markers, marker_coords = _load_virtual_markers()
marker_offsets = _load_marker_offsets()

# add some additional "markers" based on anatomical/functional landmarks
_add_synthetic_markers(markers, marker_coords)

def get_equiv_vmarker_coords(fw_name):
    """
    Return the coordinates of the gait2392 virtual marker coordinate
    equivalent to the fieldwork landmark name
    """
    g2392_name = marker_name_map[fw_name]
    return marker_coords[g2392_name]
