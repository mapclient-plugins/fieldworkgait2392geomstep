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
fieldworkmodeldict : dict
    Dictionary of model names : fieldwork models of the lower limb bones

Outputs
-------
opensimmodel : opensim.model instance
    The customised gait2392 opensim model