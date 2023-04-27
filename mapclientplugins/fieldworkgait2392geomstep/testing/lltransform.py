import numpy as np


def _trim_angle(a):
    if a < -np.pi:
        return a + 2*np.pi
    elif a > np.pi:
        return a - 2*np.pi
    else:
        return a


class LLTransformData(object):
    SHAPEMODESMAX = 100

    def __init__(self):
        self._pelvisRigid = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self._hipRot = np.array([0.0, 0.0, 0.0])
        self._kneeRot = np.array([0.0, 0.0, 0.0])
        self.nShapeModes = 1
        self.shapeModes = [0, ]
        self._shapeModeWeights = np.zeros(self.SHAPEMODESMAX, dtype=float)
        self.uniformScaling = 1.0
        self.pelvisScaling = 1.0
        self.femurScaling = 1.0
        self.patellaScaling = 1.0
        self.tibfibScaling = 1.0
        self.kneeDOF = False
        self.kneeCorr = False
        self.lastTransformSet = None

        self._shapeModelX = None
        self._uniformScalingX = None
        self._perBoneScalingX = None

    @property
    def pelvis_rigid(self):
        return self._pelvisRigid

    @pelvis_rigid.setter
    def pelvis_rigid(self, value):
        if len(value) != 6:
            raise ValueError('input pelvis_rigid vector not of length 6')
        else:
            self._pelvisRigid = np.array([value[0], value[1], value[2],
                                          _trim_angle(value[3]),
                                          _trim_angle(value[4]),
                                          _trim_angle(value[5]),
                                          ])

    @property
    def hip_rot(self):
        return self._hipRot

    @hip_rot.setter
    def hip_rot(self, value):
        if len(value) != 3:
            raise ValueError('input hip_rot vector not of length 3')
        else:
            self._hipRot = np.array([_trim_angle(v) for v in value])

    @property
    def knee_rot(self):
        if self.kneeDOF:
            return self._kneeRot[[0, 2]]
        else:
            return self._kneeRot[[0]]

    @knee_rot.setter
    def knee_rot(self, value):
        if self.kneeDOF:
            self._kneeRot[0] = _trim_angle(value[0])
            self._kneeRot[2] = _trim_angle(value[1])
        else:
            self._kneeRot[0] = _trim_angle(value[0])
    
    @property
    def shape_mode_weights(self):
        return self._shapeModeWeights[:self.nShapeModes]

    @shape_mode_weights.setter
    def shape_mode_weights(self, value):
        self._shapeModeWeights[:len(value)] = value

    # gets a flat array, sets using a list of arrays.
    @property
    def shape_model_x(self):
        self._shapeModelX = np.hstack([
                                self.shape_mode_weights[:self.nShapeModes],
                                self.pelvis_rigid,
                                self.hip_rot,
                                self.knee_rot
                                ])
        return self._shapeModelX

    @shape_model_x.setter
    def shape_model_x(self, value):
        self._shapeModelX = value
        self.shape_mode_weights = value[0]
        self.pelvis_rigid = value[1]
        self.hip_rot = value[2]
        self.knee_rot = value[3]
        self.lastTransformSet = self.shape_model_x

    @property
    def uniform_scaling_x(self):
        self._uniformScalingX = np.hstack([
                                self.uniformScaling,
                                self.pelvis_rigid,
                                self.hip_rot,
                                self.knee_rot
                                ])
        return self._uniformScalingX

    @uniform_scaling_x.setter
    def uniform_scaling_x(self, value):
        print(value)
        self._uniformScalingX = value
        self.uniformScaling = value[0]
        self.pelvis_rigid = value[1]
        self.hip_rot = value[2]
        self.knee_rot = value[3]

        # propagate isotropic scaling to each bone
        self.pelvisScaling = self.uniformScaling
        self.femurScaling = self.uniformScaling
        self.patellaScaling = self.uniformScaling
        self.tibfibScaling = self.uniformScaling

        self.lastTransformSet = self.uniform_scaling_x

    @property
    def per_bone_scaling_x(self):
        self._perBoneScalingX = np.hstack([
                                self.pelvisScaling,
                                self.femurScaling,
                                self.patellaScaling,
                                self.tibfibScaling,
                                self.pelvis_rigid,
                                self.hip_rot,
                                self.knee_rot
                                ])
        return self._perBoneScalingX

    @per_bone_scaling_x.setter
    def per_bone_scaling_x(self, value):
        self._perBoneScalingX = value
        self.pelvisScaling = value[0][1][0]
        self.femurScaling = value[0][1][1]
        self.patellaScaling = value[0][1][2]
        self.tibfibScaling = value[0][1][3]
        self.pelvis_rigid = value[1]
        self.hip_rot = value[2]
        self.knee_rot = value[3]
        self.lastTransformSet = self.per_bone_scaling_x
