"""
Module for reading Gait2392's virtual marker set.
"""

import os
import opensim
import numpy as np

from gias3.musculoskeletal.osim import Marker

SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
MARKERSET_PATH = str(os.path.join(SELF_DIR, 'data',
                                  'gait2392_Scale_MarkerSet.xml'))
MARKER_OFFSET_PATH = str(os.path.join(SELF_DIR, 'data/', 'marker_offsets.dat'))

# Dictionary mapping fieldwork landmark names to gait2392 virtual marker names
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
    virtual_markers = {}
    virtual_marker_coords = {}

    _osim_markerset = opensim.MarkerSet(MARKERSET_PATH)
    for mi in range(_osim_markerset.getSize()):
        osim_marker = _osim_markerset.get(mi)
        marker = Marker(osim_marker.clone())
        virtual_markers[marker.name] = marker
        virtual_marker_coords[marker.name] = marker.location

    return virtual_markers, virtual_marker_coords


def _add_synthetic_markers(markers_in, marker_coords_in):
    # Add hip joint centres to markers (joint coords taken from gait2392).
    markers_in['L.HJC'] = Marker(name='L.HJC', frame_name='pelvis',
                                 location=(-0.0707, -0.0661, -0.0835))
    marker_coords_in['L.HJC'] = markers_in['L.HJC'].location
    markers_in['R.HJC'] = Marker(name='R.HJC', frame_name='pelvis',
                                 location=(-0.0707, -0.0661, 0.0835))
    marker_coords_in['R.HJC'] = markers_in['R.HJC'].location

    # Add femur head centre to markers.
    markers_in['L.FHC'] = Marker(name='L.FHC', frame_name='femur_l',
                                 location=(0, 0, 0))
    marker_coords_in['L.FHC'] = markers_in['L.FHC'].location
    markers_in['R.FHC'] = Marker(name='R.FHC', frame_name='femur_r',
                                 location=(0, 0, 0))
    marker_coords_in['R.FHC'] = markers_in['R.FHC'].location

    # Add knee centre in tibia frame to markers.
    markers_in['L.Tib.KJC'] = Marker(name='L.Tib.KJC', frame_name='tibia_l',
                                     location=(0, 0, 0))
    marker_coords_in['L.Tib.KJC'] = markers_in['L.Tib.KJC'].location
    markers_in['R.Tib.KJC'] = Marker(name='R.Tib.KJC', frame_name='tibia_r',
                                     location=(0, 0, 0))
    marker_coords_in['R.Tib.KJC'] = markers_in['R.Tib.KJC'].location


def _load_marker_offsets():
    """
    Read from file the offsets between the opensim virtual markers and 
    their corresponding fieldwork landmarks
    """
    marker_offsets_list = {}
    with open(MARKER_OFFSET_PATH, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line[0] == '#':
            pass
        else:
            words = line.split()
            if len(words) == 4:
                marker_offsets_list[words[0]] = np.array(
                    [float(x) for x in words[1:]])

    return marker_offsets_list


# Load the virtual markers from file.
markers, marker_coords = _load_virtual_markers()
marker_offsets = _load_marker_offsets()
g2392_markers = markers.copy()

# Add some additional "markers" based on anatomical/functional landmarks.
_add_synthetic_markers(markers, marker_coords)


def get_equiv_vmarker_coords(fw_name):
    """
    Return the coordinates of the gait2392 virtual marker coordinate equivalent
    to the fieldwork landmark name
    """
    g2392_name = marker_name_map[fw_name]
    return marker_coords[g2392_name]
