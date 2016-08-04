
import os
from PySide import QtGui
from mapclientplugins.fieldworkgait2392geomstep.ui_configuredialog import Ui_Dialog
from mapclientplugins.fieldworkgait2392geomstep.gait2392geomcustomiser import VALID_UNITS, VALID_MODEL_MARKERS
from mapclientplugins.fieldworkgait2392geomstep.landmarktablewidget import LandmarkComboBoxTextTable


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

        # table of model and input marker pairs
        self.markerTable = LandmarkComboBoxTextTable(
                                VALID_MODEL_MARKERS,
                                self._ui.tableWidgetLandmarks,
                                )

        self._setupDialog()
        self._makeConnections()

    def _setupDialog(self):
        for s in VALID_UNITS:
            self._ui.comboBox_in_unit.addItem(s)
            self._ui.comboBox_out_unit.addItem(s)

        self._ui.lineEdit_subject_mass.setValidator(QtGui.QDoubleValidator())

    def _makeConnections(self):
        self._ui.lineEdit_id.textChanged.connect(self.validate)
        self._ui.lineEdit_osim_output_dir.textChanged.connect(self._osimOutputDirEdited)
        self._ui.pushButton_osim_output_dir.clicked.connect(self._osimOutputDirClicked)
        self._ui.pushButton_addLandmark.clicked.connect(self.markerTable.addLandmark)
        self._ui.pushButton_removeLandmark.clicked.connect(self.markerTable.removeLandmark)

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
        if osimOutputDirValid:
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
        config['in_unit'] = self._ui.comboBox_in_unit.currentText()
        config['out_unit'] = self._ui.comboBox_out_unit.currentText()
        config['adj_marker_pairs'] = self.markerTable.getLandmarkPairs()
        print('DING')
        print(config['adj_marker_pairs'])
        
        subject_mass = str(self._ui.lineEdit_subject_mass.text())
        if len(subject_mass)==0 or (subject_mass is None):
            config['subject_mass'] = None
        else:
            config['subject_mass'] = float(subject_mass)

        if self._ui.checkBox_preserve_mass_dist.isChecked():
            config['preserve_mass_distribution'] = True
        else:
            config['preserve_mass_distribution'] = False

        if self._ui.checkBox_write_osim_file.isChecked():
            config['write_osim_file'] = True
        else:
            config['write_osim_file'] = False

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

        self._ui.comboBox_in_unit.setCurrentIndex(
            VALID_UNITS.index(
                config['in_unit']
                )
            )
        self._ui.comboBox_out_unit.setCurrentIndex(
            VALID_UNITS.index(
                config['out_unit']
                )
            )

        for mm, im in sorted(config['adj_marker_pairs'].items()):
            self.markerTable.addLandmark(mm, im)

        if config['subject_mass'] is not None:
            self._ui.lineEdit_subject_mass.setText(str(config['subject_mass']))

        if config['preserve_mass_distribution']:
            self._ui.checkBox_preserve_mass_dist.setChecked(bool(True))
        else:
            self._ui.checkBox_preserve_mass_dist.setChecked(bool(False))

        if config['write_osim_file']:
            self._ui.checkBox_write_osim_file.setChecked(bool(True))
        else:
            self._ui.checkBox_write_osim_file.setChecked(bool(False))

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
        if location:
            self._previousOsimOutputDir = location
            self._ui.lineEdit_osim_output_dir.setText(location)

    def _osimOutputDirEdited(self):
        self.validate()
