import os
from PySide2 import QtGui, QtWidgets
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

        self._previousOsimOutputDir = ''

        # Table of model and input marker pairs
        self.markerTable = LandmarkComboBoxTextTable(VALID_MODEL_MARKERS, self._ui.tableWidgetLandmarks)

        self._setupDialog()
        self._makeConnections()

    def _setupDialog(self):
        for s in VALID_UNITS:
            self._ui.comboBox_in_unit.addItem(s)
            self._ui.comboBox_out_unit.addItem(s)

        self._ui.lineEdit_subject_mass.setValidator(QtGui.QDoubleValidator())

    def _makeConnections(self):
        self._ui.lineEdit_osim_output_dir.textChanged.connect(self._osimOutputDirEdited)
        self._ui.pushButton_osim_output_dir.clicked.connect(self._osimOutputDirClicked)
        self._ui.pushButton_addLandmark.clicked.connect(self.markerTable.addLandmark)
        self._ui.pushButton_removeLandmark.clicked.connect(self.markerTable.removeLandmark)

    def setWorkflowLocation(self, location):
        self._workflow_location = location

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
        location_valid = os.path.exists(os.path.join(self._workflow_location, self._ui.lineEdit_osim_output_dir.text()))
        self._ui.lineEdit_osim_output_dir.setStyleSheet(DEFAULT_STYLE_SHEET if location_valid else INVALID_STYLE_SHEET)

        self._ui.buttonBox.button(QtWidgets.QDialogButtonBox.Ok).setEnabled(location_valid)

        return location_valid

    def getConfig(self):
        """
        Get the current value of the configuration from the dialog.
        """
        config = {}
        config['osim_output_dir'] = self._ui.lineEdit_osim_output_dir.text()
        config['in_unit'] = self._ui.comboBox_in_unit.currentText()
        config['out_unit'] = self._ui.comboBox_out_unit.currentText()
        config['adj_marker_pairs'] = self.markerTable.getLandmarkPairs()

        print('DING')
        print(config['adj_marker_pairs'])

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
        Set the current value of the configuration for the dialog.
        """
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
        location = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select Directory', self._previousOsimOutputDir)
        if location:
            self._previousOsimOutputDir = location
            self._ui.lineEdit_osim_output_dir.setText(os.path.relpath(
                location, self._workflow_location))

    def _osimOutputDirEdited(self):
        self.validate()
