Fieldwork Gait2392 Geom Step
================================
MAP Client plugin for customising the OpenSim Gait2392 model using
a GIAS2 lowerlimb model (see fieldworklowerlimb2sidegenerationstep).

The following opensim model components are modified to match the input
lowerlimb model.

- Bodies
    - Mass, inertia, centre of mass : scaled according to the difference
        between the reference Gait2392 markers and markers on the input model.
    - Mass can be normalised to given subject mass or simply scaled to given subject mass while preserving reference mass distribution.
    - Display mesh for OpenSim visualtion are generated from the input model.
- Joints
    - Location and LocationInParent according to the position of joints in the input model.
    - Knee joint splines are customised to match the tibia trajectory with respect to knee flexion in the input model.
- Muscles
    - Path points are **not** modified. See fieldworkgait2392musclehmfstep for a plugin for modifying path points.
    - Wrapping objects are **not** supported.
    - Tendon slack lengths and optimal fibre lengths are updated to reflect new joint positions.
    
Body and Joint parameters of bodies not in the GIAS2 lowerlimb model (torso, talus, calcaneus, toes) are scale by the average scale factor of the model bodies. See the Scale Factor Calculation section more detail about how scale factors are calculated.

Requires
--------
- GIAS2: https://bitbucket.org/jangle/gias2
- MAP Client: https://github.com/MusculoskeletalAtlasProject/mapclient

Inputs
------
- **gias-lowerlimb** [GIAS2 LowerLimbAtlas instance]: GIAS2 lowerlimb model to be used to customise Gait2392.
- **fieldworkmodeldict** [dict][optional] :
    Fieldwork bone models to be used to customise Gait2392. The models will be combined in a GIAS2 LowerLimbAtlas model.
    Dictionary keys should be
        pelvis,
        femur-l,
        femur-r,
        patella-l,
        patella-r,
        tibiafibula-l,
        and tibiafibula-r.

Outputs
-------
- **opensimmodel** [opensim.Model instance] : The customised Gait2392 opensim model.
- **gias-lowerlimb** [GIAS2 LowerLimbAtlas instance]: Lowerlimb model used in the customisation.

Configuration
-------------
- **identifier** : Unique name for the step.
- **Output Directory** : Path of directory to output modified opensim model files.
- **Input Unit** : Unit of measurement in the input lowerlimb model or fieldwork models
- **Output Unit** : Unit of measurement of the customised OpenSim model. Default is in metres.
- **Output Osim File** : Whether to write out the customised OpenSim model.
- **Scale Other Bodies** : Whether to scale bodies not in the input lowerlimb model (torso, talus, calcaneus, toes).
- **Subject Mass** : Mass of the subject (optional). If a value is given, the scaled masses of bodies will be normalised so that total mass equals the subject mass.
- **Preserve Mass Distribution** : Whether to preserve the mass distribution amongst bodies in the reference Gait2392 model. If checked, subject mass must be given, and body masses will be simply scaled by the ratio subject_mass/ref_mass where ref_mass is the total mass of bodies in the reference Gait2392 model.

Scale Factor Calculation
------------------------
Orthotropic scale factors are calculated for each Gait2392 body to scale body mass and inertia parameters.
In addition, for bodies without a corresponding GIAS2 LowerlimbAtlas model bone (torso, talus, calcaneus, toes), the scale factors are also used to scale joint locations.

In general, given the distance l_d between 2 default Gait2392 markers and the distance l_s between the 2 corresponding markers on the GIAS2 lowerlimb model, the scale factors is l_s/l_d.
An offset is applied to the model landmarks to account for the difference between the reference markers which are on the skin surface, and input model markers which are on the bone surface.
The offset for each marker was calculated by customising the GIAS2 lowerlimb model to the reference Gait2392 model and calculating the vectors from each GIAS2 model landmark to the corresponding Gait2392 marker.

Different markers are used to calculate the scale factors in difference bodies and in different anatomical directions.
For each of the GIAS2 lowerlimb model bodies, their scale factors are calculated thus:

- Pelvis
    - x scaling: Sacral to midpoint(LASIS, RASIS) distance
    - y scaling: average of x and z scaling
    - z scaling: LASIS to RASIS distance
- Femur
    - x scaling: average of y and z scaling
    - y scaling: centre of femoral head to midpoint(MEC, LEC) distance
    - z scaling: MEC to LEC distance
- Tibia
    - x scaling: average of y and z scaling
    - y scaling: midpoint(LEC,MEC) to midpoint(LM, MM) distance
    - z scaling: MM to LM distance
    
Torso, talus, calcaneus, and toe scaling factors are the average of the scaling factors for the three bodies above.


