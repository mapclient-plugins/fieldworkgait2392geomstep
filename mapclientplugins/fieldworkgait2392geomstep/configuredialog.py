import os
from PySide6 import QtGui, QtWidgets
from mapclientplugins.fieldworkgait2392geomstep.ui_configuredialog import Ui_Dialog
from mapclientplugins.fieldworkgait2392geomstep.gait2392geomcustomiser import VALID_UNITS, VALID_MODEL_MARKERS
from mapclientplugins.fieldworkgait2392geomstep.landmarktablewidget import LandmarkComboBoxTextTable

INVALID_STYLE_SHEET = 'background-color: rgba(239, 0, 0, 50)'
DEFAULT_STYLE_SHEET = ''


class ConfigureDialog(QtWidgets.QDialog):
    """
    Configure dialog to present the user with the options to configure this step.
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        QtWidgets.QDialog.__init__(self, parent)

        self._ui = Ui_Dialog()
        self._ui.setupUi(self)

        self._workflow_location = None

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
        self._ui.lineEdit_osim_output_dir.textChanged.connect(
            self._osimOutputDirEdited)
        self._ui.pushButton_osim_output_dir.clicked.connect(
            self._osimOutputDirClicked)
        self._ui.pushButton_addLandmark.clicked.connect(
            self.markerTable.addLandmark)
        self._ui.pushButton_removeLandmark.clicked.connect(
            self.markerTable.removeLandmark)

    def setWorkflowLocation(self, location):
        self._workflow_location = location

    def _output_location(self, location=None):
        if location is None:
            display_path = self._ui.lineEdit_osim_output_dir.text()
        else:
            display_path = location
        if self._workflow_location and os.path.isabs(display_path):
            display_path = os.path.relpath(display_path, self._workflow_location)

        return display_path

    def accept(self):
        """
        Override the accept method so that we can confirm saving an
        invalid configuration.
        """
        result = QtWidgets.QMessageBox.Yes
        if not self.validate():
            result = QtWidgets.QMessageBox.warning(
                self, 'Invalid Configuration',
                'This configuration is invalid. Unpredictable behaviour may '
                'result if you choose \'Yes\', are you sure you want to save '
                'this configuration?)',
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No)

        if result == QtWidgets.QMessageBox.Yes:
            QtWidgets.QDialog.accept(self)

    def validate(self):
        """
        Validate the configuration dialog fields.  For any field that is not
        valid set the style sheet to the INVALID_STYLE_SHEET.  Return the
        outcome of the overall validity of the configuration.
        """
        # Determine if the current identifier is unique throughout the workflow. The identifierOccursCount method is part of the interface
        # to the workflow framework.
        id_value = self.identifierOccursCount(self._ui.lineEdit_id.text())
        id_valid = (id_value == 0) or (id_value == 1 and self._previousIdentifier == self._ui.lineEdit_id.text())
        self._ui.lineEdit_id.setStyleSheet(DEFAULT_STYLE_SHEET if id_valid else INVALID_STYLE_SHEET)

        file_path = self._output_location()
        if self._workflow_location:
            file_path = os.path.join(self._workflow_location, file_path)

        location_valid = os.path.exists(file_path)
        self._ui.lineEdit_osim_output_dir.setStyleSheet(DEFAULT_STYLE_SHEET if location_valid else INVALID_STYLE_SHEET)

        valid = id_valid and location_valid
        self._ui.buttonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(valid)

        return valid

    def getConfig(self):
        """
        Get the current value of the configuration from the dialog.  Also
        set the _previousIdentifier value so that we can check uniqueness of
        the identifier over the whole of the workflow.
        """
        self._previousIdentifier = self._ui.lineEdit_id.text()
        config = {}
        config['identifier'] = self._ui.lineEdit_id.text()
        config['osim_output_dir'] = self._ui.lineEdit_osim_output_dir.text()
        config['in_unit'] = self._ui.comboBox_in_unit.currentText()
        config['out_unit'] = self._ui.comboBox_out_unit.currentText()
        config['adj_marker_pairs'] = self.markerTable.getLandmarkPairs()
        # print('DING')
        # print(config['adj_marker_pairs'])

        subject_mass = str(self._ui.lineEdit_subject_mass.text())
        if len(subject_mass) == 0 or (subject_mass is None):
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
        """
        Set the current value of the configuration for the dialog.  Also
        set the _previousIdentifier value so that we can check uniqueness of
        the identifier over the whole of the workflow.
        """
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
        location = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Directory', self._previousOsimOutputDir)

        if location:
            self._previousOsimOutputDir = location
            display_location = self._output_location(location)
            self._ui.lineEdit_osim_output_dir.setText(display_location)

    def _osimOutputDirEdited(self):
        self.validate()
