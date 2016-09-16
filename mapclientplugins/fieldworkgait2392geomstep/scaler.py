"""
Module to calculate scale factors for segments in the Gait2392 model.

For each segment, (x,y,z) scale factors are calculated from differences in
distances between landmarks/markers on the customised geometry and default
2392 model.
"""

import numpy as np
from numpy.linalg import inv
import virtualmarker
from gias2.musculoskeletal import osim
from gias2.common import math

def _dist(x1, x2):
    return np.sqrt(((x1-x2)**2.0).sum(-1))

def _apply_marker_offset(model, name, coords):
    """
    Apply an offset to landmark coordinates. The offset is from the bone
    surface landmark to the skin surface opensim virtual marker.
    """

    offset_local = virtualmarker.marker_offsets[name]
    # if no offset
    if np.all(offset_local==np.zeros(3)):
        return coords

    offset_mag = math.mag(offset_local)
    offset_v = math.norm(offset_local)
    offset_v_global = math.norm(
        np.dot(
            inv(model.acs.local_transform)[:3,:3],
            offset_v[:,np.newaxis]
            ).squeeze()
        )
    offset_global = offset_v_global*offset_mag
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

#========#
# Pelvis #
#========#
# landmarks used: LASIS, RASIS, Sacral
# x scaling: Sacral to midpoint(LASIS, RASIS) distance
# y scaling: average of x and z scaling (?)
# z scaling: LASIS to RASIS distance
def calc_pelvis_scale_factors(ll, unitscale, offset=True):

    # get customised model landmarks
    cust_LASIS = _get_cust_landmark(ll.models['pelvis'], 'pelvis-LASIS', offset)
    cust_RASIS = _get_cust_landmark(ll.models['pelvis'], 'pelvis-RASIS', offset)
    cust_sacral = _get_cust_landmark(ll.models['pelvis'], 'pelvis-Sacral', offset)
    cust_lhjc = _get_cust_landmark(ll.models['pelvis'], 'pelvis-LHJC', offset)
    cust_rhjc = _get_cust_landmark(ll.models['pelvis'], 'pelvis-RHJC', offset)
    cust_o = 0.5*(cust_LASIS + cust_RASIS)
    cust_ydist = 0.5*(_dist(cust_lhjc, cust_LASIS)+_dist(cust_rhjc, cust_RASIS))

    # get reference model landmarks
    ref_LASIS = virtualmarker.get_equiv_vmarker_coords('pelvis-LASIS')
    ref_RASIS = virtualmarker.get_equiv_vmarker_coords('pelvis-RASIS')
    ref_sacral = virtualmarker.get_equiv_vmarker_coords('pelvis-Sacral')
    ref_lhjc = virtualmarker.get_equiv_vmarker_coords('pelvis-LHJC')
    ref_rhjc = virtualmarker.get_equiv_vmarker_coords('pelvis-RHJC')
    ref_o = 0.5*(ref_LASIS + ref_RASIS)
    ref_ydist = 0.5*(_dist(ref_lhjc, ref_LASIS)+_dist(ref_rhjc, ref_RASIS))

    # calculate scaling factors
    sf_x = unitscale*_dist(cust_sacral, cust_o)/_dist(ref_sacral, ref_o)
    sf_y = unitscale*cust_ydist/ref_ydist
    sf_z = unitscale*_dist(cust_LASIS, cust_RASIS)/_dist(ref_LASIS, ref_RASIS)
    # sf_y = 0.5*(sf_x + sf_z)
    
    # print('pelvis scaling factor: {:5.2f} {:5.2f} {:5.2f}'.format(sf_x, sf_y, sf_z))

    return np.array([sf_x, sf_y, sf_z])


#=======#
# Femur #
#=======#
# landmarks used: LEC, MEC, HC
# x scaling: average of y and z (?)
# y scaling: head to midpoint(MEC, LEC) distance
# z scaling: MEC to LEC distance
def calc_femur_scale_factors(ll, unitscale, side=None, offset=True):
    if side is None:
        sf_l = _calc_femur_scale_factors(ll, unitscale, 'l', offset)
        sf_r = _calc_femur_scale_factors(ll, unitscale, 'r', offset)
        return (sf_l+sf_r)*0.5
    else:
        return _calc_femur_scale_factors(ll, unitscale, side, offset)

def _calc_femur_scale_factors(ll, unitscale, side='l', offset=True):

    # get customised model landmarks
    cust_LEC = _get_cust_landmark(ll.models['femur-'+side], 'femur-LEC-'+side, offset)
    cust_MEC = _get_cust_landmark(ll.models['femur-'+side], 'femur-MEC-'+side, offset)
    cust_HC = _get_cust_landmark(ll.models['femur-'+side], 'femur-HC-'+side, offset)
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

    # print('femur scaling factor: {:5.2f} {:5.2f} {:5.2f}'.format(sf_x, sf_y, sf_z))

    return np.array([sf_x, sf_y, sf_z])

#=======#
# Tibia #
#=======#
# landmarks used: LM, MM, (LEC, MEC/KJC in tibia frame)
# x scaling: average of y and z (?)
# y scaling: midpoint(LEC,MEC) to midpoint(LM, MM) distance
# z scaling: MM to LM distance
def calc_tibia_scale_factors(ll, unitscale, side=None, offset=True):
    if side is None:
        sf_l = _calc_tibia_scale_factors(ll, unitscale, 'l', offset)
        sf_r = _calc_tibia_scale_factors(ll, unitscale, 'r', offset)
        return (sf_l+sf_r)*0.5
    else:
        return _calc_tibia_scale_factors(ll, unitscale, side, offset)

def _calc_tibia_scale_factors(ll, unitscale, side='l', offset=True):

    # get customised model landmarks
    cust_LM = _get_cust_landmark(ll.models['tibiafibula-'+side], 'tibiafibula-LM-'+side, offset)
    cust_MM = _get_cust_landmark(ll.models['tibiafibula-'+side], 'tibiafibula-MM-'+side, offset)
    cust_LEC = _get_cust_landmark(ll.models['femur-'+side], 'femur-LEC-'+side, offset)
    cust_MEC = _get_cust_landmark(ll.models['femur-'+side], 'femur-MEC-'+side, offset)
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

    # print('tibia scaling factor: {:5.2f} {:5.2f} {:5.2f}'.format(sf_x, sf_y, sf_z))

    return np.array([sf_x, sf_y, sf_z])

#============#
# whole body #
#============#
# average of scaling factors from the 3 bodies above
# to be applied to non-atlas segments e.g. torso, feet
def calc_whole_body_scale_factors(ll, unitscale, offset=True):

    sf_pelvis = calc_pelvis_scale_factors(ll, unitscale, offset)
    sf_femur_l = calc_femur_scale_factors(ll, unitscale, 'l', offset)
    sf_femur_r = calc_femur_scale_factors(ll, unitscale, 'r', offset)
    sf_tibia_l = calc_tibia_scale_factors(ll, unitscale, 'l', offset)
    sf_tibia_r = calc_tibia_scale_factors(ll, unitscale, 'r', offset)
    sf_all = np.array([
        sf_pelvis,
        sf_femur_l,
        sf_femur_r,
        sf_tibia_l,
        sf_tibia_r,
        ])

    av_sf = sf_all.mean(0)
    # print('body scaling factor: {:5.2f} {:5.2f} {:5.2f}'.format(*av_sf))
    return av_sf

#==========================#
# Scaling helper functions #
#==========================#
def _get_segment_muscles(model, segname):
    """
    Return osim muscles instances of muscles in the defined segment
    """
    seg_muscles = []
    for mus in model.muscles.values():
        # check each path point to see if they are on the segment
        for pp in mus.getAllPathPoints():
            if pp.body.name==segname:
                seg_muscles.append(mus)
                break

    return seg_muscles

def _get_segments_muscles(model, segnames):
    """
    Return osim muscles instances of muscles in the defined segments
    """
    seg_muscles = []
    for mus in model.muscles.values():
        # check each path point to see if they are on the segment
        for pp in mus.getAllPathPoints():
            if pp.body.name in segnames:
                seg_muscles.append(mus)
                break

    return seg_muscles

def scale_body_mass_inertia(body, sf):
    """
    Scales a body's inertial and mass properties according to a
    3-tuple of scale factors. Uses body scaling to calculate the
    new parameters but reverts scaling to 1 and sets parameters 
    explicitly.

    inputs
    ------
    body: osim.Body instance
    sf: 3-tuple of x,y,z scale factors

    returns
    -------
    body: scaled osim.Body instance
    """
    # print('scaling {} by {}'.format(body.name, sf))

    old_mass = body.mass
    old_mass_center = body.massCenter
    old_inertia = body.inertia

    # apply opensim scaling
    # body.scale(sf, False)
    body.scaleInertialProperties(sf, True)

    # get scaled mass and inertial tensor
    new_mass = body.mass
    new_mass_center = body.massCenter
    new_inertia = body.inertia

    # scale back to 1,1,1
    sf_inv = 1.0/np.array(sf)
    # body.scale(sf_inv, False)
    body.scaleInertialProperties(sf_inv, True)

    # set new params
    body.mass = new_mass
    body.massCenter = new_mass_center
    body.inertia = new_inertia

    print('mass: {} -> {}'.format(old_mass, new_mass))
    print('mass_center: {} -> {}'.format(old_mass_center, new_mass_center))
    print('inertia: {} -> {}'.format(old_inertia.diagonal(), new_inertia.diagonal()))

    return body

def scale_joint(joint, sfs, bodies):
    """
    Scales a joints properties according to a 3-tuple of scale factors.
    Uses joint.scale to calculate the new parameters but reverts scaling
    to 1 and sets parameters explicitly.

    inputs
    ------
    joint: osim.Joint instance
    sfs: a list of 3-tuple of x,y,z scale factors
    bodies: list of names of bodies connected to joint

    returns
    -------
    joint: scaled osim.Joint instance
    """
    print('scaling {} by {}'.format(joint.name, sfs))

    old_location = joint.location
    old_locationInParent = joint.locationInParent

    # apply opensim scaling
    joint_scales = []
    for bi, bname in enumerate(bodies):
        s = osim.Scale(
                sfs[bi],
                '{}_scale_{}'.format(joint.name, bname),
                bname,
                )
        joint_scales.append(s)

    joint.scale(*joint_scales)

    # get new joint properties
    new_location = joint.location
    new_locationInParent = joint.locationInParent

    # scale back to 1,1,1
    # sf_inv = 1.0/np.array(sf)
    joint_scales = []
    for bi, bname in enumerate(bodies):
        s = osim.Scale(
                1.0/sfs[bi],
                '{}_scale_{}'.format(joint.name, bname),
                bname,
                )
        joint_scales.append(s)

    joint.scale(*joint_scales)

    # set new params
    joint.location = new_location
    joint.locationInParent = new_locationInParent

    print('location: {} -> {}'.format(old_location, new_location))
    print('locationInParent: {} -> {}'.format(old_locationInParent, new_locationInParent))

    return joint

def scale_body_muscle_params(model, bodyscales, state):
    """
    Scales a body's muscle tendon slack length and optimal fibre length
    properties according to a 3-tuple of scale factors. Uses Muscle.scale()
    to calculate the new parameters but reverts scaling to 1 and 
    sets parameters explicitly.

    inputs
    ------
    model: osim.Model instance
    bodyscales: dictionary of body names and 3-tuple body scale factors
    state: model state at which scaling is to be performed

    returns
    -------
    None
    """

    raise(DeprecationWarning, 'This does not work as it does not use muscle prescale, scale and postscale')

    # get all muscles of the body
    muscles = _get_segments_muscles(model, list(bodyscales.keys()))
    if len(muscles)==0:
        print('WARNING: no muscles found for bodies {}'.format(list(bodyscales.keys())))

    # for each muscle, apply scaling, get opt fibre length and
    # tendon slack length, scale back, then apply those two params
    for mus in muscles:
        print('Scaling muscle {}'.format(mus.name))
        old_ofl = mus.optimalFiberLength
        old_tsl = mus.tendonSlackLength

        mus_scales = []
        for bname, bscale in bodyscales.items():
            s = osim.Scale(
                    bscale,
                    '{}_scale_{}'.format(mus.name, bname),
                    bname,
                    )
            mus_scales.append(s)
        mus.scale(state, *mus_scales)

        new_ofl = mus.optimalFiberLength
        new_tsl = mus.tendonSlackLength


        mus_inv_scales = []
        for bname, bscale in bodyscales.items():
            s = osim.Scale(
                    1.0/np.array(bscale),
                    '{}_inv_scale_{}'.format(mus.name, bname),
                    bname,
                    )
            mus_inv_scales.append(s)
        mus.scale(state, *mus_inv_scales)

        mus.optimalFiberLength = new_ofl
        mus.tendonSlackLength = new_tsl

        print('optimal fiber length: {} -> {}'.format(old_ofl, new_ofl))
        print('tendon slack length: {} -> {}'.format(old_tsl, new_tsl))

def calc_scale_factors_all_bodies(LL, unit_scaling, scale_other_bodies=True):
    """
    Returns a list of scale factors, on for each body in the model
    """

    sf_list = []

    # pelvis
    sf_list.append(
        osim.Scale(
            calc_pelvis_scale_factors(
                LL, unit_scaling,
                ),
            'pelvis_scaling', 'pelvis'
            )
        )

    # femur
    sf_list.append(
        osim.Scale(
            calc_femur_scale_factors(
                LL, unit_scaling, None,
                ),
            'femur_l_scaling', 'femur_l'
            )
        )
    sf_list.append(
        osim.Scale(
            calc_femur_scale_factors(
                LL, unit_scaling, None,
                ),
            'femur_r_scaling', 'femur_r'
            )
        )

    # tibia
    sf_list.append(
        osim.Scale(
            calc_tibia_scale_factors(
                LL, unit_scaling, None,
                ),
            'tibia_l_scaling', 'tibia_l'
            )
        )
    sf_list.append(
        osim.Scale(
            calc_tibia_scale_factors(
                LL, unit_scaling, None,
                ),
            'tibia_r_scaling', 'tibia_r'
            )
        )

    if scale_other_bodies:
        sf_whole = calc_whole_body_scale_factors(
            LL, unit_scaling
            )

        # torso
        sf_list.append(
            osim.Scale(sf_whole, 'torso_scaling', 'torso')
            )

        # talus
        sf_list.append(
            osim.Scale(sf_whole, 'talus_l_scaling', 'talus_l')
            )
        sf_list.append(
            osim.Scale(sf_whole, 'talus_r_scaling', 'talus_r')
            )

        # calcn
        sf_list.append(
            osim.Scale(sf_whole, 'calcn_l_scaling', 'calcn_l')
            )
        sf_list.append(
            osim.Scale(sf_whole, 'calcn_r_scaling', 'calcn_r')
            )

        # toes
        sf_list.append(
            osim.Scale(sf_whole, 'toes_l_scaling', 'toes_l')
            )
        sf_list.append(
            osim.Scale(sf_whole, 'toes_r_scaling', 'toes_r')
            )

    return sf_list

def scale_body_visual_geometry(body, sf):

    _sf = osim.opensim.Vec3(sf[0], sf[1], sf[2])
    gset = body._osimBody.getDisplayer().getGeometrySet()
    for gi in range(gset.getSize()):
        g = gset.get(gi)
        g.setScaleFactors(_sf)

