"""
Module for reading Gait2392's virtual markerset
"""

import os
import opensim
import numpy as np

SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
MARKERSET_PATH = os.path.join(SELF_DIR, 'data', 'gait2392_Scale_MarkerSet.xml')

# dictionary mapping fieldwork landmark names to gait2392's virtual marker
# names
marker_name_map = {
    'pelvis-RASIS': 'R.ASIS',
    'pelvis-LASIS': 'L.ASIS',
    'pelvis-Sacral': 'V.Sacral',
    'femur-LEC-r': 'R.Knee.Lat',
    'femur-MEC-r': 'R.Knee.Med',
    'tibiafibula-LM-r': 'R.Ankle.Lat',
    'tibiafibula-MM-r': 'R.Ankle.Mat',
    'femur-LEC-l': 'L.Knee.Lat',
    'femur-MEC-l': 'L.Knee.Med',
    'tibiafibula-LM-l': 'L.Ankle.Lat',
    'tibiafibula-MM-l': 'L.Ankle.Mat',
    'femur-HC-l': 'L.HJC',
    'femur-HC-l': 'R.HJC',
    'tibiafibula-KJC-l': 'L.Tib.KJC',
    'tibiafibula-KJC-r': 'R.Tib.KJC',
}

class Marker(object):
    """
    Pythonic wrap of opensim's Marker class
    """

    def __init__(self, m=None, bodyname=None, offset=None):
        if m is None:
            self._osimMarker = opensim.Marker()
            self.bodyName = bodyname
            self.offset = offset
        else:
            self._osimMarker = m

    @property
    def name(self):
        return self._osimMarker.getName()

    @name.setter
    def name(self, name):
        self._osimMarker.setName(name)

    @property
    def bodyName(self):
        return self._osimMarker.getBodyName()

    @bodyName.setter
    def bodyName(self, bodyName):
        self._osimMarker.setBodyName(bodyName)

    @property
    def offset(self):
        v = opensim.Vec3()
        self._osimMarker.getOffset(v)
        return np.array([v.get(i) for i in range(3)])

    @offset.setter
    def offset(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimMarker.setOffset(v)

    # Same as location
    @property
    def offset(self):
        v = opensim.Vec3()
        self._osimMarker.getOffset(v)
        return np.array([v.get(i) for i in range(3)])

    # Same as location
    @offset.setter
    def offset(self, x):
        v = opensim.Vec3(x[0], x[1], x[2])
        self._osimMarker.setOffset(v)

def _load_virtual_markers():
    markers = {}
    marker_coords = {}
    _osim_markerset = opensim.MarkerSet(MARKERSET_PATH)
    for mi in range(_osim_markerset.getSize()):
        osim_marker = _osim_markerset.get(mi)
        marker = Marker(osim_marker)
        markers[marker.name] = marker
        marker_coords[marker.name] = marker.offset

    return markers, marker_coords

# load up vitual markers from file
markers, marker_coords = _load_virtual_markers()

# add femur head centre to markers
markers['L.HJC'] = Marker(bodyname='femur_l', offset=(0,0,0))
marker_coords['L.HJC'] = markers['L.HJC'].offset
markers['R.HJC'] = Marker(bodyname='femur_r', offset=(0,0,0))
marker_coords['R.HJC'] = markers['R.HJC'].offset

# add knee centre in tibia frame to markers
markers['L.Tib.KJC'] = Marker(bodyname='tibia_l', offset=(0,0,0))
marker_coords['L.Tib.KJC'] = markers['L.Tib.KJC'].offset
markers['R.Tib.KJC'] = Marker(bodyname='tibia_r', offset=(0,0,0))
marker_coords['R.Tib.KJC'] = markers['R.Tib.KJC'].offset

def get_equiv_vmarker_coords(fw_name):
    g2392_name = marker_name_map[fw_name]
    return marker_coords[g2392_name]