import numpy as np

def _trimAngle(a):
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
        self.shapeModes = [0,]
        self._shapeModeWeights = np.zeros(self.SHAPEMODESMAX, dtype=float)
        self.uniformScaling = 1.0
        self.pelvisScaling = 1.0
        self.femurScaling = 1.0
        self.petallaScaling = 1.0
        self.tibfibScaling = 1.0
        self.kneeDOF = False
        self.kneeCorr = False
        self.lastTransformSet = None

        self._shapeModelX = None
        self._uniformScalingX = None
        self._perBoneScalingX = None

    @property
    def pelvisRigid(self):
        return self._pelvisRigid

    @pelvisRigid.setter
    def pelvisRigid(self, value):
        if len(value)!=6:
            raise ValueError('input pelvisRigid vector not of length 6')
        else:
            self._pelvisRigid = np.array([value[0], value[1], value[2],
                                          _trimAngle(value[3]),
                                          _trimAngle(value[4]),
                                          _trimAngle(value[5]),
                                         ])

    @property
    def hipRot(self):
        return self._hipRot

    @hipRot.setter
    def hipRot(self, value):
        if len(value)!=3:
            raise ValueError('input hipRot vector not of length 3')
        else:
            self._hipRot = np.array([_trimAngle(v) for v in value])

    @property
    def kneeRot(self):
        if self.kneeDOF:
            return self._kneeRot[[0,2]]
        else:
            return self._kneeRot[[0]]

    @kneeRot.setter
    def kneeRot(self, value):
        if self.kneeDOF:
            self._kneeRot[0] = _trimAngle(value[0])
            self._kneeRot[2] = _trimAngle(value[1])
        else:
            self._kneeRot[0] = _trimAngle(value[0])
    
    @property
    def shapeModeWeights(self):
        return self._shapeModeWeights[:self.nShapeModes]

    @shapeModeWeights.setter
    def shapeModeWeights(self, value):
        self._shapeModeWeights[:len(value)] = value

    # gets a flat array, sets using a list of arrays.
    @property
    def shapeModelX(self):
        self._shapeModelX = np.hstack([
                                self.shapeModeWeights[:self.nShapeModes],
                                self.pelvisRigid,
                                self.hipRot,
                                self.kneeRot
                                ])
        return self._shapeModelX

    @shapeModelX.setter
    def shapeModelX(self, value):
        a = self.nShapeModes
        self._shapeModelX = value
        self.shapeModeWeights = value[0]
        self.pelvisRigid = value[1]
        self.hipRot = value[2]
        self.kneeRot = value[3]
        self.lastTransformSet = self.shapeModelX

    @property
    def uniformScalingX(self):
        self._uniformScalingX = np.hstack([
                                self.uniformScaling,
                                self.pelvisRigid,
                                self.hipRot,
                                self.kneeRot
                                ])
        return self._uniformScalingX

    @uniformScalingX.setter
    def uniformScalingX(self, value):
        print value
        a = 1
        self._uniformScalingX = value
        self.uniformScaling = value[0]
        self.pelvisRigid = value[1]
        self.hipRot = value[2]
        self.kneeRot = value[3]

        # propagate isotropic scaling to each bone
        self.pelvisScaling = self.uniformScaling
        self.femurScaling = self.uniformScaling
        self.patellaScaling = self.uniformScaling
        self.tibfibScaling = self.uniformScaling

        self.lastTransformSet = self.uniformScalingX

    @property
    def perBoneScalingX(self):
        self._perBoneScalingX = np.hstack([
                                self.pelvisScaling,
                                self.femurScaling,
                                self.patellaScaling,
                                self.tibfibScaling,
                                self.pelvisRigid,
                                self.hipRot,
                                self.kneeRot
                                ])
        return self._perBoneScalingX

    @perBoneScalingX.setter
    def perBoneScalingX(self, value):
        a = 4
        self._perBoneScalingX = value
        self.pelvisScaling = value[0][1][0]
        self.femurScaling = value[0][1][1]
        self.patellaScaling = value[0][1][2]
        self.tibfibScaling = value[0][1][3]
        self.pelvisRigid = value[1]
        self.hipRot = value[2]
        self.kneeRot = value[3]
        self.lastTransformSet = self.perBoneScalingX
