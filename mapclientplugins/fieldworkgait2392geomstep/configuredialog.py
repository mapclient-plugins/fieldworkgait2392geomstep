

from PySide import QtGui
from mapclientplugins.fieldworkgait2392geomstep.ui_configuredialog import Ui_Dialog

INVALID_STYLE_SHEET = 'background-color: rgba(239, 0, 0, 50)'
DEFAULT_STYLE_SHEET = ''

REG_MODES = ('shapemodel', 'uniformscaling', 'perbonescaling', 'manual')
SIDEOPTIONS = ('left', 'right')

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
        # Set a place holder for a callable that will get set from the step.
        # We will use this method to decide whether the identifier is unique.
        self.identifierOccursCount = None

        self._makeConnections()

        for regmode in REG_MODES:
            self._ui.comboBox_regmode.addItem(regmode)

        for side in SIDEOPTIONS:
            self._ui.comboBox_side.addItem(side)

    def _makeConnections(self):
        self._ui.lineEdit_id.textChanged.connect(self.validate)

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
        value = self.identifierOccursCount(self._ui.lineEdit_id.text())
        valid = (value == 0) or (value == 1 and self._previousIdentifier == self._ui.lineEdit_id.text())
        if valid:
            self._ui.lineEdit_id.setStyleSheet(DEFAULT_STYLE_SHEET)
        else:
            self._ui.lineEdit_id.setStyleSheet(INVALID_STYLE_SHEET)

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
        config['registration_mode'] = self._ui.comboBox_regmode.currentText()
        config['side'] = self._ui.comboBox_side.currentText()
        config['pcs_to_fit'] = str(self._ui.spinBox_pcsToFit.value())
        config['mweight'] = str(self._ui.doubleSpinBox_mWeight.value())
        config['pelvis-RASIS'] = self._ui.lineEdit_RASIS.text()
        config['pelvis-LASIS'] = self._ui.lineEdit_LASIS.text()
        config['pelvis-Sacral'] = self._ui.lineEdit_Sacral.text()
        config['femur-LEC'] = self._ui.lineEdit_LEC.text()
        config['femur-MEC'] = self._ui.lineEdit_MEC.text()
        config['tibiafibula-LM'] = self._ui.lineEdit_LM.text()
        config['tibiafibula-MM'] = self._ui.lineEdit_MM.text()
        config['marker_radius'] = self._ui.doubleSpinBox_markerRadius.value()
        config['skin_pad'] = self._ui.doubleSpinBox_skinPad.value()
        if self._ui.checkBox_kneecorr.isChecked():
            config['knee_corr'] = 'True'
        else:
            config['knee_corr'] = 'False'
        if self._ui.checkBox_kneedof.isChecked():
            config['knee_dof'] = 'True'
        else:
            config['knee_dof'] = 'False'
        if self._ui.checkBox_GUI.isChecked():
            config['GUI'] = 'True'
        else:
            config['GUI'] = 'False'
        return config

    def setConfig(self, config):
        '''
        Set the current value of the configuration for the dialog.  Also
        set the _previousIdentifier value so that we can check uniqueness of the
        identifier over the whole of the workflow.
        '''
        self._previousIdentifier = config['identifier']
        self._ui.lineEdit_id.setText(config['identifier'])
        self._ui.comboBox_regmode.setCurrentIndex(
            REG_MODES.index(config['registration_mode'])
            )
        if config.get('side'):
            self._ui.comboBox_side.setCurrentIndex(
                SIDEOPTIONS.index(config['side'])
                )
        else:
            self._ui.comboBox_side.setCurrentIndex(0)

        self._ui.spinBox_pcsToFit.setValue(int(config['pcs_to_fit']))
        self._ui.doubleSpinBox_mWeight.setValue(float(config['mweight']))
        self._ui.lineEdit_RASIS.setText(config['pelvis-RASIS'])
        self._ui.lineEdit_LASIS.setText(config['pelvis-LASIS'])
        self._ui.lineEdit_Sacral.setText(config['pelvis-Sacral'])
        self._ui.lineEdit_LEC.setText(config['femur-LEC'])
        self._ui.lineEdit_MEC.setText(config['femur-MEC'])
        self._ui.lineEdit_MM.setText(config['tibiafibula-MM'])
        self._ui.lineEdit_LM.setText(config['tibiafibula-LM'])
        self._ui.doubleSpinBox_markerRadius.setValue(float(config['marker_radius']))
        self._ui.doubleSpinBox_skinPad.setValue(float(config['skin_pad']))
        if config['knee_corr']=='True':
            self._ui.checkBox_kneecorr.setChecked(bool(True))
        else:
            self._ui.checkBox_kneecorr.setChecked(bool(False))
        if config['knee_dof']=='True':
            self._ui.checkBox_kneedof.setChecked(bool(True))
        else:
            self._ui.checkBox_kneedof.setChecked(bool(False))
        if config['GUI']=='True':
            self._ui.checkBox_GUI.setChecked(bool(True))
        else:
            self._ui.checkBox_GUI.setChecked(bool(False))

