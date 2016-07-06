"""
Module to calculate scale factors for segments in the Gait2392 model.

For each segment, (x,y,z) scale factors are calculated from differences in
distances between landmarks/markers on the customised geometry and default
2392 model.
"""

import numpy as np
import virtualmarker

def _dist(x1, x2):
    return np.sqrt(((x1-x2)**2.0).sum(-1))

#========#
# Pelvis #
#========#
# landmarks used: LASIS, RASIS, Sacral
# x scaling: Sacral to midpoint(LASIS, RASIS) distance
# y scaling: average of x and z scaling (?)
# z scaling: LASIS to RASIS distance
def calc_pelvis_scale_factors(ll, unitscale):

    # get customised model landmarks
    cust_LASIS = ll.models['pelvis'].landmarks['pelvis-LASIS']
    cust_RASIS = ll.models['pelvis'].landmarks['pelvis-RASIS']
    cust_sacral = ll.models['pelvis'].landmarks['pelvis-Sacral']
    cust_o = 0.5*(cust_LASIS + cust_RASIS)

    # get reference model landmarks
    ref_LASIS = virtualmarker.get_equiv_vmarker_coords('pelvis-LASIS')
    ref_RASIS = virtualmarker.get_equiv_vmarker_coords('pelvis-RASIS')
    ref_sacral = virtualmarker.get_equiv_vmarker_coords('pelvis-Sacral')
    ref_o = 0.5*(ref_LASIS + ref_RASIS)

    # calculate scaling factors
    sf_x = unitscale*_dist(cust_sacral, cust_o)/_dist(ref_sacral, ref_o)
    sf_z = unitscale*_dist(cust_LASIS, cust_RASIS)/_dist(ref_LASIS, ref_RASIS)
    sf_y = 0.5*(sf_x + sf_z)

    return np.array([sf_x, sf_y, sf_z])


#=======#
# Femur #
#=======#
# landmarks used: LEC, MEC, HC
# x scaling: average of y and z (?)
# y scaling: head to midpoint(MEC, LEC) distance
# z scaling: MEC to LEC distance
def calc_femur_scale_factors(ll, unitscale, side='l'):

    # get customised model landmarks
    cust_LEC = ll.models['femur-'+side].landmarks['femur-LEC']
    cust_MEC = ll.models['femur-'+side].landmarks['femur-MEC']
    cust_HC = ll.models['femur-'+side].landmarks['femur-HC']
    cust_o = 0.5*(cust_LEC + cust_MEC)

    # get reference model landmarks
    ref_LEC = virtualmarker.get_equiv_vmarker_coords('femur-LEC-'+side)
    ref_MEC = virtualmarker.get_equiv_vmarker_coords('femur-MEC-'+side)
    ref_HC = virtualmarker.get_equiv_vmarker_coords('femur-HC-'+side)
    ref_o = 0.5*(ref_LEC + ref_MEC)

    # calculate scaling factors
    sf_y = unitscale*_dist(cust_HC, cust_o)/_dist(ref_HC, ref_o)
    sf_z = unitscale*_dist(cust_LEC, cust_MEC)/_dist(ref_LEC, ref_MEC)
    sf_x = 0.5*(sf_y + sf_z)

    return np.array([sf_x, sf_y, sf_z])

#=======#
# Tibia #
#=======#
# landmarks used: LM, MM, (LEC, MEC/KJC in tibia frame)
# x scaling: average of y and z (?)
# y scaling: midpoint(LEC,MEC) to midpoint(LM, MM) distance
# z scaling: MM to LM distance
def calc_tibia_scale_factors(ll, unitscale, side='l'):

    # get customised model landmarks
    cust_LM = ll.models['tibiafibula-'+side].landmarks['tibiafibula-LM']
    cust_MM = ll.models['tibiafibula-'+side].landmarks['tibiafibula-MM']
    cust_LEC = ll.models['femur-'+side].landmarks['femur-LEC']
    cust_MEC = ll.models['femur-'+side].landmarks['femur-MEC']
    cust_kjc = 0.5*(cust_LEC + cust_MEC)
    cust_ajc = 0.5*(cust_LM + cust_MM)

    # get reference model landmarks
    ref_LM = virtualmarker.get_equiv_vmarker_coords('tibiafibula-LM-'+side)
    ref_MM = virtualmarker.get_equiv_vmarker_coords('tibiafibula-MM-'+side)
    ref_kjc = virtualmarker.get_equiv_vmarker_coords('tibiafibula-KJC-'+side)
    ref_ajc = 0.5*(ref_LM + ref_MM)

    # calculate scaling factors
    sf_y = unitscale*_dist(cust_kjc, cust_ajc)/_dist(ref_kjc, ref_ajc)
    sf_z = unitscale*_dist(cust_LM, cust_MM)/_dist(ref_LM, ref_MM)
    sf_x = 0.5*(sf_y + sf_z)

    return np.array([sf_x, sf_y, sf_z])

#============#
# whole body #
#============#
# average of scaling factors from the 3 bodies above
# to be applied to non-atlas segments e.g. torso, feet
def calc_whole_body_scale_factors(ll, unitscale):

    sf_pelvis = calc_pelvis_scale_factors(ll, unitscale)
    sf_femur_l = calc_femur_scale_factors(ll, unitscale, 'l')
    sf_femur_r = calc_femur_scale_factors(ll, unitscale, 'r')
    sf_tibia_l = calc_tibia_scale_factors(ll, unitscale, 'l')
    sf_tibia_r = calc_tibia_scale_factors(ll, unitscale, 'r')
    sf_all = np.array([
        sf_pelvis,
        sf_femur_l,
        sf_femur_r,
        sf_tibia_l,
        sf_tibia_r,
        ])

    av_sf = sf_all.mean(0)
    return av_sf