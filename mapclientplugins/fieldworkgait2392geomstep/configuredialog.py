
import os
from PySide import QtGui
from mapclientplugins.fieldworkgait2392geomstep.ui_configuredialog import Ui_Dialog

INVALID_STYLE_SHEET = 'background-color: rgba(239, 0, 0, 50)'
DEFAULT_STYLE_SHEET = ''

class ConfigureDialog(QtGui.QDialog):
    '''
    Configure dialog to present the user with the options to configure this step.
    '''

    def __init__(self, parent=None):
        '''
        Constructor
        '''
        QtGui.QDialog.__init__(self, parent)
        
        self._ui = Ui_Dialog()
        self._ui.setupUi(self)

        # Keep track of the previous identifier so that we can track changes
        # and know how many occurrences of the current identifier there should
        # be.
        self._previousIdentifier = ''
        self._previousOsimOutputDir = ''
        # Set a place holder for a callable that will get set from the step.
        # We will use this method to decide whether the identifier is unique.
        self.identifierOccursCount = None

        self._makeConnections()

    def _makeConnections(self):
        self._ui.lineEdit_id.textChanged.connect(self.validate)
        self._ui.lineEdit_osim_output_dir.textChanged.connect(self._osimOutputDirEdited)
        self._ui.pushButton_osim_output_dir.clicked.connect(self._osimOutputDirClicked)

    def accept(self):
        '''
        Override the accept method so that we can confirm saving an
        invalid configuration.
        '''
        result = QtGui.QMessageBox.Yes
        if not self.validate():
            result = QtGui.QMessageBox.warning(self, 'Invalid Configuration',
                'This configuration is invalid.  Unpredictable behaviour may result if you choose \'Yes\', are you sure you want to save this configuration?)',
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if result == QtGui.QMessageBox.Yes:
            QtGui.QDialog.accept(self)

    def validate(self):
        '''
        Validate the configuration dialog fields.  For any field that is not valid
        set the style sheet to the INVALID_STYLE_SHEET.  Return the outcome of the 
        overall validity of the configuration.
        '''
        # Determine if the current identifier is unique throughout the workflow
        # The identifierOccursCount method is part of the interface to the workflow framework.
        idValue = self.identifierOccursCount(self._ui.lineEdit_id.text())
        idValid = (idValue == 0) or (idValue == 1 and self._previousIdentifier == self._ui.lineEdit_id.text())
        if idValid:
            self._ui.lineEdit_id.setStyleSheet(DEFAULT_STYLE_SHEET)
        else:
            self._ui.lineEdit_id.setStyleSheet(INVALID_STYLE_SHEET)

        osimOutputDirValid = os.path.exists(self._ui.lineEdit_osim_output_dir.text())
        if fileLocValid:
            self._ui.lineEdit_osim_output_dir.setStyleSheet(DEFAULT_STYLE_SHEET)
        else:
            self._ui.lineEdit_osim_output_dir.setStyleSheet(INVALID_STYLE_SHEET)
            
        valid = idValid and osimOutputDirValid
        self._ui.buttonBox.button(QtGui.QDialogButtonBox.Ok).setEnabled(valid)

        return valid

    def getConfig(self):
        '''
        Get the current value of the configuration from the dialog.  Also
        set the _previousIdentifier value so that we can check uniqueness of the
        identifier over the whole of the workflow.
        '''
        self._previousIdentifier = self._ui.lineEdit_id.text()
        config = {}
        config['identifier'] = self._ui.lineEdit_id.text()
        config['osim_output_dir'] = self._ui.lineEdit_osim_output_dir.text()
        
        if self._ui.checkBox_write_osim_file.isChecked():
            config['write_osim_file'] = True
        else:
            config['write_osim_file'] = False

        if self._ui.checkBox_convert_mm_to_m.isChecked():
            config['convert_mm_to_m'] = True
        else:
            config['convert_mm_to_m'] = False

        if self._ui.checkBox_scale_other_bodies.isChecked():
            config['scale_other_bodies'] = True
        else:
            config['scale_other_bodies'] = False

        if self._ui.checkBox_GUI.isChecked():
            config['GUI'] = True
        else:
            config['GUI'] = False

        return config

    def setConfig(self, config):
        '''
        Set the current value of the configuration for the dialog.  Also
        set the _previousIdentifier value so that we can check uniqueness of the
        identifier over the whole of the workflow.
        '''
        self._previousIdentifier = config['identifier']
        self._ui.lineEdit_id.setText(config['identifier'])
        self._previousOsimOutputDir = config['osim_output_dir']
        self._ui.lineEdit_osim_output_dir.setText(config['osim_output_dir'])

        if config['write_osim_file']:
            self._ui.checkBox_write_osim_file.setChecked(bool(True))
        else:
            self._ui.checkBox_write_osim_file.setChecked(bool(False))

        if config['convert_mm_to_m']:
            self._ui.checkBox_convert_mm_to_m.setChecked(bool(True))
        else:
            self._ui.checkBox_convert_mm_to_m.setChecked(bool(False))

        if config['scale_other_bodies']:
            self._ui.checkBox_scale_other_bodies.setChecked(bool(True))
        else:
            self._ui.checkBox_scale_other_bodies.setChecked(bool(False))

        if config['GUI']:
            self._ui.checkBox_GUI.setChecked(bool(True))
        else:
            self._ui.checkBox_GUI.setChecked(bool(False))

    def _osimOutputDirClicked(self):
        location = QtGui.QFileDialog.getExistingDirectory(self, 'Select Directory', self._previousOsimOutputDir)
        if location[0]:
            self._previousOsimOutputDir = location[0]
            self._ui.lineEdit_osim_output_dir.setText(location[0])

    def _osimOutputDirEdited(self):
        self.validate()
