Fieldwork Gait2392 Geom Step
================================
MAP Client plugin for customising the OpenSim Gait2392 model geometry using
fieldwork models. Parameters modified are body frame definitions, visual
meshes, and the scaling of non-patient-specific bodies.

Requires
--------
GIAS2: https://bitbucket.org/jangle/gias2,
VTK, Mayavi2

Inputs
------
gias-lowerlimb : gias2.musculoskeletal.bonemodel.LowerLimbAtlas instance
    Lower limb model to be used to customise gait2392.
fieldworkmodeldict : dict [optional]
    Bone models to be used to customisation gait2392.
    Dictionary keys should be:
        pelvis
        femur-l
        femur-r
        patella-l
        patella-r
        tibiafibula-l
        tibiafibula-r

Outputs
-------
opensimmodel : opensim.model instance
    The customised gait2392 opensim model
gias-lowerlimb : gias2.musculoskeletal.bonemodel.LowerLimbAtlas instance
    The lowerlimb model used in the customisation
