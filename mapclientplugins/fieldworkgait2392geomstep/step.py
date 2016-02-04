
'''
MAP Client Plugin Step
'''
import os
import json

from PySide import QtGui

from mapclient.mountpoints.workflowstep import WorkflowStepMountPoint
from mapclientplugins.fieldworkgait2392geomstep.configuredialog import ConfigureDialog

from mapclientplugins.fieldworkgait2392geomstep import llstep
from mapclientplugins.fieldworkgait2392geomstep.fieldworkgait2392dialog import FieldworkGait2392GeomDialog

SELF_DIR = os.path.split(os.path.realpath(__file__))[0]
TEMPLATE_OSIM_FILENAME = os.path.join(SELF_DIR, 'data/gait2392_simbody.osim')

class FieldworkGait2392GeomStep(WorkflowStepMountPoint):
    '''
    Step for customising the OpenSim Gait2392 model geometry using
    fieldwork models. Parameters modified are body frame definitions, visual
    meshes, and the scaling of non-patient-specific bodies.

    Inputs
    ------
    fieldworkmodeldict : dict
        Dictionary of model names : fieldwork models of the lower limb bones

    Outputs
    -------
    opensimmodel : opensim.model instance
        The customised gait2392 opensim model
    '''

    def __init__(self, location):
        super(FieldworkGait2392GeomStep, self).__init__('Fieldwork Gait2392 Geometry Customisation', location)
        self._configured = False # A step cannot be executed until it has been configured.
        self._category = 'OpenSim'
        # Add any other initialisation code here:
        # Ports:
        self.addPort(('http://physiomeproject.org/workflow/1.0/rdf-schema#port',
                      'http://physiomeproject.org/workflow/1.0/rdf-schema#uses',
                      'ju#fieldworkmodeldict'))
        self.addPort(('http://physiomeproject.org/workflow/1.0/rdf-schema#port',
                      'http://physiomeproject.org/workflow/1.0/rdf-schema#provides',
                      'http://physiomeproject.org/workflow/1.0/rdf-schema#osimmodel'))

        self._config = {}
        self._config['identifier'] = ''
        self._config['GUI'] = 'True'
        self._config['scale_other_bodies'] = False
        self._config
        for l in LLLANDMARKS:
            self._config[l] = ''

        self._data = llstep.LLStepData(self._config)

    def execute(self):
        '''
        Add your code here that will kick off the execution of the step.
        Make sure you call the _doneExecution() method when finished.  This method
        may be connected up to a button in a widget for example.
        '''
        # Put your execute step code here before calling the '_doneExecution' method.
        self._data.loadData()
        self._data.updateFromConfig()
        print('LL estimation configs:')
        print self._data.config
        if self._config['GUI']=='True':
            # start gui
            self._widget = LowerLimbGenerationDialog(self._data, self._doneExecution)
            self._widget.setModal(True)
            self._setCurrentWidget(self._widget)
        else:
            self._data.register()
            self._doneExecution()

    def setPortData(self, index, dataIn):
        '''
        Add your code here that will set the appropriate objects for this step.
        The index is the index of the port in the port list.  If there is only one
        uses port for this step then the index can be ignored.
        '''
        if index == 0:
            self._data.inputLandmarks = dataIn # http://physiomeproject.org/workflow/1.0/rdf-schema#landmarks
        elif index == 1:
            self._data.inputPCs = dataIn # ju#principalcomponents
        else:
            self._data.inputModelDict = dataIn # ju#fieldworkmodeldict

    def getPortData(self, index):
        '''
        Add your code here that will return the appropriate objects for this step.
        The index is the index of the port in the port list.  If there is only one
        provides port for this step then the index can be ignored.
        '''
        if index == 3:
            return self._data.outputModelDict
        else:
            return self._data.outputTransform

    def configure(self):
        '''
        This function will be called when the configure icon on the step is
        clicked.  It is appropriate to display a configuration dialog at this
        time.  If the conditions for the configuration of this step are complete
        then set:
            self._configured = True
        '''
        dlg = ConfigureDialog()
        dlg.identifierOccursCount = self._identifierOccursCount
        dlg.setConfig(self._config)
        dlg.validate()
        dlg.setModal(True)
        
        if dlg.exec_():
            self._config = dlg.getConfig()
            self._data.config = self._config
        
        self._configured = dlg.validate()
        self._configuredObserver()

    def getIdentifier(self):
        '''
        The identifier is a string that must be unique within a workflow.
        '''
        return self._config['identifier']

    def setIdentifier(self, identifier):
        '''
        The framework will set the identifier for this step when it is loaded.
        '''
        self._config['identifier'] = identifier

    def serialize(self):
        '''
        Add code to serialize this step to disk. Returns a json string for
        mapclient to serialise.
        '''
        return json.dumps(self._config, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def deserialize(self, string):
        '''
        Add code to deserialize this step from disk. Parses a json string
        given by mapclient
        '''
        self._config.update(json.loads(string))

        d = ConfigureDialog()
        d.identifierOccursCount = self._identifierOccursCount
        d.setConfig(self._config)
        self._configured = d.validate()


