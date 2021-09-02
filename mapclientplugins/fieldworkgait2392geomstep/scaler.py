"""
Module to calculate scale factors for segments in the Gait2392 model.

For each segment, (x,y,z) scale factors are calculated from differences in
distances between landmarks/markers on the customised geometry and default
2392 model.
"""

import numpy as np
from numpy.linalg import inv
from gias2.musculoskeletal import osim
from gias2.common import math

# from mapclientplugins.fieldworkgait2392geomstep import virtualmarker
from . import virtualmarker


def _dist(x1, x2):
    return np.sqrt(((x1 - x2) ** 2.0).sum(-1))


def _apply_marker_offset(model, name, coords):
    """
    Apply an offset to landmark coordinates. The offset is from the bone
    surface landmark to the skin surface opensim virtual marker.
    """

    offset_local = virtualmarker.marker_offsets[name]
    # if no offset
    if np.all(offset_local == np.zeros(3)):
        return coords

    offset_mag = math.mag(offset_local)
    offset_v = math.norm(offset_local)
    offset_v_global = math.norm(
        np.dot(
            inv(model.acs.local_transform)[:3, :3],
            offset_v[:, np.newaxis]
        ).squeeze()
    )
    offset_global = offset_v_global * offset_mag
    return coords + offset_global


def _get_cust_landmark(body, lname, offset=True):
    if lname[-2:] in ('-l', '-r'):
        _lname = lname[:-2]
    else:
        _lname = lname
    ld = body.landmarks[_lname]
    if offset:
        # apply an offset from bone surface to skin surface
        return _apply_marker_offset(body, lname, ld)
    else:
        return ld


# ========#
# Pelvis #
# ========#
# landmarks used: LASIS, RASIS, Sacral
# x scaling: Sacral to midpoint(LASIS, RASIS) distance
# y scaling: average of x and z scaling (?)
# z scaling: LASIS to RASIS distance
def calc_pelvis_scale_factors(ll, unitscale, offset=True):
    # Get customised model landmarks.
    cust_lasis = _get_cust_landmark(ll.models['pelvis'], 'pelvis-LASIS',
                                    offset)
    cust_rasis = _get_cust_landmark(ll.models['pelvis'], 'pelvis-RASIS',
                                    offset)
    cust_sacral = _get_cust_landmark(ll.models['pelvis'], 'pelvis-Sacral',
                                     offset)
    cust_lhjc = _get_cust_landmark(ll.models['pelvis'], 'pelvis-LHJC', offset)
    cust_rhjc = _get_cust_landmark(ll.models['pelvis'], 'pelvis-RHJC', offset)
    cust_o = 0.5 * (cust_lasis + cust_rasis)
    cust_ydist = 0.5 * (_dist(cust_lhjc, cust_lasis) + _dist(cust_rhjc,
                                                             cust_rasis))

    # Get reference model landmarks.
    ref_lasis = virtualmarker.get_equiv_vmarker_coords('pelvis-LASIS')
    ref_rasis = virtualmarker.get_equiv_vmarker_coords('pelvis-RASIS')
    ref_sacral = virtualmarker.get_equiv_vmarker_coords('pelvis-Sacral')
    ref_lhjc = virtualmarker.get_equiv_vmarker_coords('pelvis-LHJC')
    ref_rhjc = virtualmarker.get_equiv_vmarker_coords('pelvis-RHJC')
    ref_o = 0.5 * (ref_lasis + ref_rasis)
    ref_ydist = 0.5 * (_dist(ref_lhjc, ref_lasis) + _dist(ref_rhjc, ref_rasis))

    # Calculate scaling factors.
    sf_x = unitscale * _dist(cust_sacral, cust_o) / _dist(ref_sacral, ref_o)
    sf_y = unitscale * cust_ydist / ref_ydist
    sf_z = unitscale * _dist(cust_lasis, cust_rasis) / _dist(ref_lasis,
                                                             ref_rasis)

    return np.array([sf_x, sf_y, sf_z])


# =======#
# Femur #
# =======#
# landmarks used: LEC, MEC, HC
# x scaling: average of y and z (?)
# y scaling: head to midpoint(MEC, LEC) distance
# z scaling: MEC to LEC distance
def calc_femur_scale_factors(ll, unitscale, side=None, offset=True):
    if side is None:
        sf_l = _calc_femur_scale_factors(ll, unitscale, 'l', offset)
        sf_r = _calc_femur_scale_factors(ll, unitscale, 'r', offset)
        return (sf_l + sf_r) * 0.5
    else:
        return _calc_femur_scale_factors(ll, unitscale, side, offset)


def _calc_femur_scale_factors(ll, unitscale, side='l', offset=True):
    # Get customised model landmarks.
    cust_lec = _get_cust_landmark(ll.models['femur-' + side],
                                  'femur-LEC-' + side, offset)
    cust_mec = _get_cust_landmark(ll.models['femur-' + side],
                                  'femur-MEC-' + side, offset)
    cust_hc = _get_cust_landmark(ll.models['femur-' + side],
                                 'femur-HC-' + side, offset)
    cust_o = 0.5 * (cust_lec + cust_mec)

    # Get reference model landmarks.
    ref_lec = virtualmarker.get_equiv_vmarker_coords('femur-LEC-' + side)
    ref_mec = virtualmarker.get_equiv_vmarker_coords('femur-MEC-' + side)
    ref_hc = virtualmarker.get_equiv_vmarker_coords('femur-HC-' + side)
    ref_o = 0.5 * (ref_lec + ref_mec)

    # Calculate scaling factors.
    sf_y = unitscale * _dist(cust_hc, cust_o) / _dist(ref_hc, ref_o)
    sf_z = unitscale * _dist(cust_lec, cust_mec) / _dist(ref_lec, ref_mec)
    sf_x = 0.5 * (sf_y + sf_z)

    return np.array([sf_x, sf_y, sf_z])


# =======#
# Tibia #
# =======#
# landmarks used: LM, MM, (LEC, MEC/KJC in tibia frame)
# x scaling: average of y and z (?)
# y scaling: midpoint(LEC,MEC) to midpoint(LM, MM) distance
# z scaling: MM to LM distance
def calc_tibia_scale_factors(ll, unitscale, side=None, offset=True):
    if side is None:
        sf_l = _calc_tibia_scale_factors(ll, unitscale, 'l', offset)
        sf_r = _calc_tibia_scale_factors(ll, unitscale, 'r', offset)
        return (sf_l + sf_r) * 0.5
    else:
        return _calc_tibia_scale_factors(ll, unitscale, side, offset)


def _calc_tibia_scale_factors(ll, unitscale, side='l', offset=True):
    # Get customised model landmarks.
    cust_lm = _get_cust_landmark(ll.models['tibiafibula-' + side],
                                 'tibiafibula-LM-' + side, offset)
    cust_mm = _get_cust_landmark(ll.models['tibiafibula-' + side],
                                 'tibiafibula-MM-' + side, offset)
    cust_lec = _get_cust_landmark(ll.models['femur-' + side],
                                  'femur-LEC-' + side, offset)
    cust_mec = _get_cust_landmark(ll.models['femur-' + side],
                                  'femur-MEC-' + side, offset)
    cust_kjc = 0.5 * (cust_lec + cust_mec)
    cust_ajc = 0.5 * (cust_lm + cust_mm)

    # Get reference model landmarks.
    ref_lm = virtualmarker.get_equiv_vmarker_coords('tibiafibula-LM-' + side)
    ref_mm = virtualmarker.get_equiv_vmarker_coords('tibiafibula-MM-' + side)
    ref_kjc = virtualmarker.get_equiv_vmarker_coords('tibiafibula-KJC-' + side)
    ref_ajc = 0.5 * (ref_lm + ref_mm)

    # Calculate scaling factors.
    sf_y = unitscale * _dist(cust_kjc, cust_ajc) / _dist(ref_kjc, ref_ajc)
    sf_z = unitscale * _dist(cust_lm, cust_mm) / _dist(ref_lm, ref_mm)
    sf_x = 0.5 * (sf_y + sf_z)

    return np.array([sf_x, sf_y, sf_z])


def calc_scale_factors_all_bodies(ll, unit_scaling, scale_other_bodies=True):
    """
    Builds and returns a list of OpenSim::Scale objects, one for each Body in
    the Model. Scale objects describe the scale_factors to be used for scaling
    the Model.

    inputs
    ------
    ll:
    unit_scaling:
    scale_other_bodies:

    outputs
    -------
    scale_list:
    """

    pelvis_scale_factors = calc_pelvis_scale_factors(ll, unit_scaling,)
    femur_l_scale_factors = calc_femur_scale_factors(ll, unit_scaling, "l",)
    femur_r_scale_factors = calc_femur_scale_factors(ll, unit_scaling, "r",)
    tibia_l_scale_factors = calc_tibia_scale_factors(ll, unit_scaling, "l",)
    tibia_r_scale_factors = calc_tibia_scale_factors(ll, unit_scaling, "r",)

    scale_list = [
        osim.Scale(pelvis_scale_factors, 'pelvis_scaling', 'pelvis'),
        osim.Scale(femur_l_scale_factors, 'femur_l_scaling', 'femur_l'),
        osim.Scale(femur_r_scale_factors, 'femur_r_scaling', 'femur_r'),
        osim.Scale(tibia_l_scale_factors, 'tibia_l_scaling', 'tibia_l'),
        osim.Scale(tibia_r_scale_factors, 'tibia_r_scaling', 'tibia_r')
    ]

    if scale_other_bodies:
        all_scale_factors = np.array([
            pelvis_scale_factors,
            femur_l_scale_factors,
            femur_r_scale_factors,
            tibia_l_scale_factors,
            tibia_r_scale_factors
        ])
        average_scale_factors = all_scale_factors.mean(0)

        # torso
        scale_list.append(
            osim.Scale(average_scale_factors, 'torso_scaling', 'torso')
        )

        # talus
        scale_list.append(
            osim.Scale(average_scale_factors, 'talus_l_scaling', 'talus_l')
        )
        scale_list.append(
            osim.Scale(average_scale_factors, 'talus_r_scaling', 'talus_r')
        )

        # calcn
        scale_list.append(
            osim.Scale(average_scale_factors, 'calcn_l_scaling', 'calcn_l')
        )
        scale_list.append(
            osim.Scale(average_scale_factors, 'calcn_r_scaling', 'calcn_r')
        )

        # toes
        scale_list.append(
            osim.Scale(average_scale_factors, 'toes_l_scaling', 'toes_l')
        )
        scale_list.append(
            osim.Scale(average_scale_factors, 'toes_r_scaling', 'toes_r')
        )

    return scale_list
