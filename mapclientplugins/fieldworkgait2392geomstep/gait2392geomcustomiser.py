"""
Gait2392 customisation
"""
import os
import numpy as np
import copy

from gias2.fieldwork.field import geometric_field
from gias2.musculoskeletal import mocap_landmark_preprocess
from gias2.musculoskeletal.bonemodels import bonemodels
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfit
from gias2.musculoskeletal.bonemodels import lowerlimbatlasfitscaling

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

SELF_DIRECTORY = os.path.split(__file__)[0]
PELVIS_SUBMESHES = ('RH', 'LH', 'sac')
PELVIS_SUBMESH_ELEMS = {'RH': range(0, 73),
                        'LH': range(73,146),
                        'sac': range(146, 260),
                        }
PELVIS_BASISTYPES = {'tri10':'simplex_L3_L3','quad44':'quad_L3_L3'}

TIBFIB_SUBMESHES = ('tibia', 'fibula')
TIBFIB_SUBMESH_ELEMS = {'tibia': range(0, 46),
                        'fibula': range(46,88),
                        }
TIBFIB_BASISTYPES = {'tri10':'simplex_L3_L3','quad44':'quad_L3_L3'}
    
class LLStepData(object):

    _shapeModelFilenameRight = os.path.join(SELF_DIRECTORY, 'data/shape_models/LLP26_right_mirrored_from_left_rigid.pc')
    _boneModelFilenamesRight = {'pelvis': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof'),
                                      os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.ens'),
                                      os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.mesh'),
                                      ),
                           'femur': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                     os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_quartic_flat.ens'),
                                     os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_right_quartic_flat.mesh'),
                                     ),
                           'patella': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_right.ens'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_right.mesh'),
                                       ),
                           'tibiafibula': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_cubic_right_mirrored_from_left_mean_rigid_LLP26.geof'),
                                           os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_right_cubic_flat.ens'),
                                           os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_right_cubic_flat.mesh'),
                                           ),
                            }
    _shapeModelFilenameLeft = os.path.join(SELF_DIRECTORY, 'data/shape_models/LLP26_rigid.pc')
    _boneModelFilenamesLeft = {'pelvis': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_mean_rigid_LLP26.geof'),
                                      os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.ens'),
                                      os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/pelvis_combined_cubic_flat.mesh'),
                                      ),
                           'femur': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_mean_rigid_LLP26.geof'),
                                     os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_quartic_flat.ens'),
                                     os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/femur_left_quartic_flat.mesh'),
                                     ),
                           'patella': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_left_mean_rigid_LLP26.geof'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_left.ens'),
                                       os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/patella_11_left.mesh'),
                                       ),
                           'tibiafibula': (os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_cubic_left_mean_rigid_LLP26.geof'),
                                           os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_left_cubic_flat.ens'),
                                           os.path.join(SELF_DIRECTORY, 'data/atlas_meshes/tibia_fibula_left_cubic_flat.mesh'),
                                           ),
                            }
    _validRegistrationModes = ('shapemodel', 'uniformscaling', 'perbonescaling')
    landmarkNames = ('pelvis-LASIS', 'pelvis-RASIS', 'pelvis-Sacral',
                      'femur-LEC', 'femur-MEC', 'tibiafibula-LM',
                      'tibiafibula-MM',
                      )
    minArgs = {'method':'L-BFGS-B',
                 'jac':False,
                 'bounds':None, 'tol':1e-6,
                 'options':{'eps':1e-5},
                 }

    def __init__(self, config):
        self.config = config
        self.T = LLTransformData()
        self.inputLandmarks = None # a dict of landmarks
        # self._targetLandmarksNames = None # list of strings matching keys in self.inputLandmarks
        self._targetLandmarks = None

        self.inputPCs = None
        self._inputModelDict = None
        self._outputModelDict = None
        self.landmarkErrors = None
        self.landmarkRMSE = None
        self.fitMDist = None

        # self.regCallback = None

    def loadData(self):
        if self.config['side']=='left':
            self.LL = bonemodels.LowerLimbLeftAtlas('lower_limb_left')
            self.LL.bone_files = self._boneModelFilenamesLeft
            self.LL.combined_pcs_filename = self._shapeModelFilenameLeft
        elif self.config['side']=='right':
            self.LL = bonemodels.LowerLimbRightAtlas('lower_limb_right')
            self.LL.bone_files = self._boneModelFilenamesRight
            self.LL.combined_pcs_filename = self._shapeModelFilenameRight

        self.LL.load_bones()

    def resetLL(self):
        self.LL.update_all_models(*self.LL._neutral_params)
        self.T = LLTransformData()
        self.landmarkErrors = None
        self.landmarkRMSE = None

    def updateFromConfig(self):
        targetLandmarkNames = [self.config[ln] for ln in self.landmarkNames]
        self.targetLandmarkNames = targetLandmarkNames
        self.nShapeModes = self.config['pcs_to_fit']
        if self.kneeCorr:
            self.LL.enable_knee_adduction_correction()
        else:
            self.LL.disable_knee_adduction_correction()
        if self.kneeDOF:
            self.LL.enable_knee_adduction_dof
        else:
            self.LL.disable_knee_adduction_dof()

    def updateLLModel(self):
        """update LL model using current transformations.
        Just shape model deformations
        """
        self.LL.update_all_models(self.T.shapeModeWeights,
                                  self.T.shapeModes,
                                  self.T.pelvisRigid,
                                  self.T.hipRot,
                                  self.T.kneeRot
                                  )

    def _preprocessLandmarks(self, l):
        return np.array(mocap_landmark_preprocess.preprocess_lower_limb(
                        self.markerRadius,
                        self.skinPad,
                        l[0], l[1], l[2], l[3], l[4], l[5], l[6]),
                        )

    @property
    def outputModelDict(self):
        self._outputModelDict = dict([(m[0], m[1].gf) for m in self.LL.models.items()])

        # add pelvis submeshes
        self._outputModelDict['pelvis flat'] = copy.deepcopy(self._outputModelDict['pelvis'])
        lh_gf, sac_gf, rh_gf = self._splitPelvisGFs()
        self._outputModelDict['hemipelvis-left'] = lh_gf
        self._outputModelDict['sacrum'] = sac_gf
        self._outputModelDict['hemipelvis-right'] = rh_gf
        # self._outputModelDict['pelvis'] = self._createNestedPelvis(self._outputModelDict['pelvis flat'])

        # add seperate tibia and fibula
        tibia_gf, fibula_gf = self._splitTibiaFibulaGFs()
        self._outputModelDict['tibia'] = tibia_gf
        self._outputModelDict['fibula'] = fibula_gf

        return self._outputModelDict

    def _splitTibiaFibulaGFs(self):
        tibfib = self.LL.models['tibiafibula'].gf
        tib = tibfib.makeGFFromElements(
                'tibia',
                TIBFIB_SUBMESH_ELEMS['tibia'],
                TIBFIB_BASISTYPES,
                )
        fib = tibfib.makeGFFromElements(
                'fibula',
                TIBFIB_SUBMESH_ELEMS['fibula'],
                TIBFIB_BASISTYPES,
                )

        return tib, fib

    def _splitPelvisGFs(self):
        """ Given a flattened pelvis model, create left hemi, sacrum,
        and right hemi meshes
        """
        gf = self.LL.models['pelvis'].gf
        lhgf = gf.makeGFFromElements(
                    'hemipelvis-left',
                    PELVIS_SUBMESH_ELEMS['LH'],
                    PELVIS_BASISTYPES
                    )
        sacgf = gf.makeGFFromElements(
                    'sacrum',
                    PELVIS_SUBMESH_ELEMS['sac'],
                    PELVIS_BASISTYPES
                    )
        rhgf = gf.makeGFFromElements(
                    'hemipelvis-right',
                    PELVIS_SUBMESH_ELEMS['RH'],
                    PELVIS_BASISTYPES
                    )
        return lhgf, sacgf, rhgf

    def _createNestedPelvis(self, gf):
        """ Given a flattened pelvis model, create a hierarchical model
        """
        newgf = geometric_field.geometric_field(
                    gf.name, 3,
                    field_dimensions=2,
                    field_basis=PELVIS_BASISTYPES
                    )

        for subname in PELVIS_SUBMESHES:
            subgf = gf.makeGFFromElements(
                        subname,
                        PELVIS_SUBMESH_ELEMS[subname],
                        PELVIS_BASISTYPES
                        )
            newgf.add_element_with_parameters(
                subgf.ensemble_field_function,
                subgf.field_parameters,
                tol=0.0
                )

        return newgf

    @property
    def outputTransform(self):
        return self.T

    @property
    def validRegistrationModes(self):
        return self._validRegistrationModes  
    
    @property
    def registrationMode(self):
        return self.config['registration_mode']
        # return self._registrationMode
    
    @registrationMode.setter
    def registrationMode(self, value):
        if value in self.validRegistrationModes:
            self.config['registration_mode'] = value
        else:
            raise ValueError('Invalid registration mode. Given {}, must be one of {}'.format(value, self.validRegistrationModes))

    @property
    def targetLandmarkNames(self):
        # return self._targetLandmarkNames
        return [self.config[ln] for ln in self.landmarkNames]

    @targetLandmarkNames.setter
    def targetLandmarkNames(self, value):
        if len(value)!=7:
            raise ValueError('7 input landmark names required for {}'.format(self._landmarkNames))
        else:
            # self._targetLandmarkNames = value
            # save to config dict
            for li, ln in enumerate(self.landmarkNames):
                self.config[ln] = value[li]
            # evaluate target landmark coordinates
            # if (self.inputLandmarks is not None) and ('' not in value):
            #     self._targetLandmarks = np.array([self.inputLandmarks[n] for n in self._targetLandmarkNames])

    @property
    def targetLandmarks(self):
        if '' in self.targetLandmarkNames:
            raise ValueError('Null string in targetLandmarkNames')

        self._targetLandmarks = np.array([self.inputLandmarks[n] for n in self.targetLandmarkNames])
        self._targetLandmarks = self._preprocessLandmarks(self._targetLandmarks)
        return self._targetLandmarks
    
    @property
    def inputModelDict(self):
        return self._inputModelDict

    @inputModelDict.setter
    def inputModelDict(self, value):
        self._inputModelDict = value

    @property
    def mWeight(self):
        return float(self.config['mweight'])

    @mWeight.setter
    def mWeight(self, value):
        self.config['mweight'] = str(value)

    @property
    def nShapeModes(self):
        return int(self.config['pcs_to_fit'])

    @property
    def markerRadius(self):
        return float(self.config['marker_radius'])

    @markerRadius.setter
    def markerRadius(self, value):
        self.config['marker_radius'] = str(value)

    @property
    def skinPad(self):
        return float(self.config['skin_pad'])

    @skinPad.setter
    def skinPad(self, value):
        self.config['skin_pad'] = str(value)
    
    @nShapeModes.setter
    def nShapeModes(self, n):
        self.config['pcs_to_fit'] = str(n)
        n = int(n)
        self.T.nShapeModes = n
        self.T.shapeModes = np.arange(n, dtype=int)

    @property
    def kneeCorr(self):
        return self.config['knee_corr']=='True'

    @kneeCorr.setter
    def kneeCorr(self, value):
        self.T.kneeCorr = value
        if value:
            self.config['knee_corr'] = 'True'
            self.LL.enable_knee_adduction_correction()
        else:
            self.config['knee_corr'] = 'False'
            self.LL.disable_knee_adduction_correction()

    @property
    def kneeDOF(self):
        return self.config['knee_dof']=='True'

    @kneeDOF.setter
    def kneeDOF(self, value):
        self.T.kneeDOF = value
        if value:
            self.config['knee_dof'] = 'True'
            self.LL.enable_knee_adduction_dof()
        else:
            self.config['knee_dof'] = 'False'
            self.LL.disable_knee_adduction_dof()

    def register(self, callbackSignal=None):
        self.updateFromConfig()
        mode = self.config['registration_mode']

        if self.targetLandmarks is None:
            raise RuntimeError('Target Landmarks not set')

        if callbackSignal is not None:
            def callback(output):
                callbackSignal.emit(output)
        else:
            callback=None

        if mode=='shapemodel':
            if self.T.shapeModes is None:
                raise RuntimeError('Number of pcs to fit not defined')
            if self.mWeight is None:
                raise RuntimeError('Mahalanobis penalty weight not defined')
            output = _registerShapeModel(self, callback)
        elif mode=='uniformscaling':
            output = _registerUniformScaling(self)
        elif mode=='perbonescaling':
            output = _registerPerBoneScaling(self)
        return output

def _registerShapeModel(lldata, callback=None):

    # if lladata.T.shapeModelX has not changed from the default,
    # use None for x0 so that it is automatically calculated
    x0Temp = lldata.T.shapeModelX
    if np.all(x0Temp==0.0):
        x0 = None
    else:
        x0 = x0Temp
    print(x0)

    # do the fit
    xFitted,\
    optLandmarkDist,\
    optLandmarkRMSE,\
    fitInfo = lowerlimbatlasfit.fit(
                    lldata.LL,
                    lldata.targetLandmarks,
                    lldata.landmarkNames,
                    lldata.T.shapeModes,
                    lldata.mWeight,
                    x0=x0,
                    minimise_args=lldata.minArgs,
                    callback=callback,
                    )
    lldata.landmarkRMSE = optLandmarkRMSE
    lldata.landmarkErrors = optLandmarkDist
    lldata.fitMDist =  fitInfo['mahalanobis_distance']
    lldata.T.shapeModelX = xFitted[-1]
    print('new X:'+str(lldata.T.shapeModelX))
    return xFitted, optLandmarkDist, optLandmarkRMSE, fitInfo

def _registerUniformScaling(lldata, callback=None):
    
    # if lladata.T.uniformScalingX has not changed from the default,
    # use None for x0 so that it is automatically calculated
    x0Temp = lldata.T.uniformScalingX
    if (x0Temp[0]==1.0) and np.all(x0Temp[1:]==0.0):
        x0 = None
    else:
        x0 = x0Temp
    print(x0)

    # do the fit
    xFitted,\
    optLandmarkDist,\
    optLandmarkRMSE,\
    fitInfo = lowerlimbatlasfitscaling.fit(
                    lldata.LL,
                    lldata.targetLandmarks,
                    lldata.landmarkNames,
                    bones_to_scale='uniform',
                    x0=x0,
                    minimise_args=lldata.minArgs,
                    # callback=callback,
                    )
    lldata.landmarkRMSE = optLandmarkRMSE
    lldata.landmarkErrors = optLandmarkDist
    lldata.fitMDist = -1.0
    lldata.T.uniformScalingX = xFitted[-1]
    print('new X:'+str(lldata.T.uniformScalingX))
    return xFitted, optLandmarkDist, optLandmarkRMSE, fitInfo

def _registerPerBoneScaling(lldata, callback=None):
    
    # if lladata.T.perboneScalingX has not changed from the default,
    # use None for x0 so that it is automatically calculated
    x0Temp = lldata.T.perBoneScalingX
    if np.all(x0Temp[:4]==1.0) and np.all(x0Temp[4:]==0.0):
        x0 = None
    else:
        x0 = x0Temp
    print(x0)
    bones = ('pelvis', 'femur', 'patella', 'tibiafibula')
    xFitted,\
    optLandmarkDist,\
    optLandmarkRMSE,\
    fitInfo = lowerlimbatlasfitscaling.fit(
                    lldata.LL,
                    lldata.targetLandmarks,
                    lldata.landmarkNames,
                    bones_to_scale=bones,
                    x0=x0,
                    minimise_args=lldata.minArgs,
                    # callback=callback,
                    )
    lldata.landmarkRMSE = optLandmarkRMSE
    lldata.landmarkErrors = optLandmarkDist
    lldata.fitMDist = -1.0
    lldata.T.perBoneScalingX = xFitted[-1]
    print('new X:'+str(lldata.T.perBoneScalingX))
    return xFitted, optLandmarkDist, optLandmarkRMSE, fitInfo