# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'configuredialog.ui'
##
## Created by: Qt User Interface Compiler version 5.15.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(639, 611)
        self.gridLayout = QGridLayout(Dialog)
        self.gridLayout.setObjectName(u"gridLayout")
        self.configGroupBox = QGroupBox(Dialog)
        self.configGroupBox.setObjectName(u"configGroupBox")
        self.formLayout = QFormLayout(self.configGroupBox)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.label = QLabel(self.configGroupBox)
        self.label.setObjectName(u"label")

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.label)

        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.lineEdit_osim_output_dir = QLineEdit(self.configGroupBox)
        self.lineEdit_osim_output_dir.setObjectName(u"lineEdit_osim_output_dir")

        self.horizontalLayout.addWidget(self.lineEdit_osim_output_dir)

        self.pushButton_osim_output_dir = QPushButton(self.configGroupBox)
        self.pushButton_osim_output_dir.setObjectName(u"pushButton_osim_output_dir")

        self.horizontalLayout.addWidget(self.pushButton_osim_output_dir)


        self.formLayout.setLayout(0, QFormLayout.FieldRole, self.horizontalLayout)

        self.label_2 = QLabel(self.configGroupBox)
        self.label_2.setObjectName(u"label_2")

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.label_2)

        self.comboBox_in_unit = QComboBox(self.configGroupBox)
        self.comboBox_in_unit.setObjectName(u"comboBox_in_unit")

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.comboBox_in_unit)

        self.label_4 = QLabel(self.configGroupBox)
        self.label_4.setObjectName(u"label_4")

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.label_4)

        self.comboBox_out_unit = QComboBox(self.configGroupBox)
        self.comboBox_out_unit.setObjectName(u"comboBox_out_unit")

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.comboBox_out_unit)

        self.label_10 = QLabel(self.configGroupBox)
        self.label_10.setObjectName(u"label_10")

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.label_10)

        self.checkBox_write_osim_file = QCheckBox(self.configGroupBox)
        self.checkBox_write_osim_file.setObjectName(u"checkBox_write_osim_file")

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.checkBox_write_osim_file)

        self.label_3 = QLabel(self.configGroupBox)
        self.label_3.setObjectName(u"label_3")

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.label_3)

        self.checkBox_scale_other_bodies = QCheckBox(self.configGroupBox)
        self.checkBox_scale_other_bodies.setObjectName(u"checkBox_scale_other_bodies")

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.checkBox_scale_other_bodies)

        self.label_5 = QLabel(self.configGroupBox)
        self.label_5.setObjectName(u"label_5")

        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.label_5)

        self.lineEdit_subject_mass = QLineEdit(self.configGroupBox)
        self.lineEdit_subject_mass.setObjectName(u"lineEdit_subject_mass")

        self.formLayout.setWidget(5, QFormLayout.FieldRole, self.lineEdit_subject_mass)

        self.label_6 = QLabel(self.configGroupBox)
        self.label_6.setObjectName(u"label_6")

        self.formLayout.setWidget(6, QFormLayout.LabelRole, self.label_6)

        self.checkBox_preserve_mass_dist = QCheckBox(self.configGroupBox)
        self.checkBox_preserve_mass_dist.setObjectName(u"checkBox_preserve_mass_dist")

        self.formLayout.setWidget(6, QFormLayout.FieldRole, self.checkBox_preserve_mass_dist)

        self.label_7 = QLabel(self.configGroupBox)
        self.label_7.setObjectName(u"label_7")

        self.formLayout.setWidget(7, QFormLayout.LabelRole, self.label_7)

        self.verticalLayout = QVBoxLayout()
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.tableWidgetLandmarks = QTableWidget(self.configGroupBox)
        if (self.tableWidgetLandmarks.columnCount() < 2):
            self.tableWidgetLandmarks.setColumnCount(2)
        __qtablewidgetitem = QTableWidgetItem()
        self.tableWidgetLandmarks.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.tableWidgetLandmarks.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        self.tableWidgetLandmarks.setObjectName(u"tableWidgetLandmarks")
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableWidgetLandmarks.sizePolicy().hasHeightForWidth())
        self.tableWidgetLandmarks.setSizePolicy(sizePolicy)
        self.tableWidgetLandmarks.horizontalHeader().setMinimumSectionSize(200)
        self.tableWidgetLandmarks.horizontalHeader().setDefaultSectionSize(200)

        self.verticalLayout.addWidget(self.tableWidgetLandmarks)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.pushButton_addLandmark = QPushButton(self.configGroupBox)
        self.pushButton_addLandmark.setObjectName(u"pushButton_addLandmark")

        self.horizontalLayout_2.addWidget(self.pushButton_addLandmark)

        self.pushButton_removeLandmark = QPushButton(self.configGroupBox)
        self.pushButton_removeLandmark.setObjectName(u"pushButton_removeLandmark")

        self.horizontalLayout_2.addWidget(self.pushButton_removeLandmark)


        self.verticalLayout.addLayout(self.horizontalLayout_2)


        self.formLayout.setLayout(7, QFormLayout.FieldRole, self.verticalLayout)

        self.label_11 = QLabel(self.configGroupBox)
        self.label_11.setObjectName(u"label_11")

        self.formLayout.setWidget(8, QFormLayout.LabelRole, self.label_11)

        self.checkBox_GUI = QCheckBox(self.configGroupBox)
        self.checkBox_GUI.setObjectName(u"checkBox_GUI")
        self.checkBox_GUI.setEnabled(False)

        self.formLayout.setWidget(8, QFormLayout.FieldRole, self.checkBox_GUI)


        self.gridLayout.addWidget(self.configGroupBox, 0, 0, 1, 1)

        self.buttonBox = QDialogButtonBox(Dialog)
        self.buttonBox.setObjectName(u"buttonBox")
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)

        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        QWidget.setTabOrder(self.lineEdit_osim_output_dir, self.pushButton_osim_output_dir)
        QWidget.setTabOrder(self.pushButton_osim_output_dir, self.comboBox_in_unit)
        QWidget.setTabOrder(self.comboBox_in_unit, self.comboBox_out_unit)
        QWidget.setTabOrder(self.comboBox_out_unit, self.checkBox_write_osim_file)
        QWidget.setTabOrder(self.checkBox_write_osim_file, self.checkBox_scale_other_bodies)
        QWidget.setTabOrder(self.checkBox_scale_other_bodies, self.lineEdit_subject_mass)
        QWidget.setTabOrder(self.lineEdit_subject_mass, self.checkBox_preserve_mass_dist)
        QWidget.setTabOrder(self.checkBox_preserve_mass_dist, self.tableWidgetLandmarks)
        QWidget.setTabOrder(self.tableWidgetLandmarks, self.pushButton_addLandmark)
        QWidget.setTabOrder(self.pushButton_addLandmark, self.pushButton_removeLandmark)
        QWidget.setTabOrder(self.pushButton_removeLandmark, self.checkBox_GUI)
        QWidget.setTabOrder(self.checkBox_GUI, self.buttonBox)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)

        QMetaObject.connectSlotsByName(Dialog)
    # setupUi

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"Configure OpenSim Gait2392 Customisation Step", None))
        self.configGroupBox.setTitle("")
        self.label.setText(QCoreApplication.translate("Dialog", u"Output Directory:", None))
        self.pushButton_osim_output_dir.setText(QCoreApplication.translate("Dialog", u"...", None))
        self.label_2.setText(QCoreApplication.translate("Dialog", u"Input Unit:", None))
        self.label_4.setText(QCoreApplication.translate("Dialog", u"Output Unit:", None))
        self.label_10.setText(QCoreApplication.translate("Dialog", u"Output .osim file:", None))
        self.checkBox_write_osim_file.setText("")
        self.label_3.setText(QCoreApplication.translate("Dialog", u"Scale other bodies:", None))
        self.checkBox_scale_other_bodies.setText("")
        self.label_5.setText(QCoreApplication.translate("Dialog", u"Subject Mass (kg):", None))
        self.label_6.setText(QCoreApplication.translate("Dialog", u"Preserve Mass Distribution:", None))
        self.checkBox_preserve_mass_dist.setText("")
        self.label_7.setText(QCoreApplication.translate("Dialog", u"Adjustable Markers:", None))
        ___qtablewidgetitem = self.tableWidgetLandmarks.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("Dialog", u"Model Markers", None));
        ___qtablewidgetitem1 = self.tableWidgetLandmarks.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("Dialog", u"Input Markers", None));
        self.pushButton_addLandmark.setText(QCoreApplication.translate("Dialog", u"Add Marker Pair", None))
        self.pushButton_removeLandmark.setText(QCoreApplication.translate("Dialog", u"Remove Marker Pair", None))
        self.label_11.setText(QCoreApplication.translate("Dialog", u"GUI:", None))
        self.checkBox_GUI.setText("")
    # retranslateUi

