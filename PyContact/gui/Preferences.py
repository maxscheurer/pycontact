# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Preferences.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_PreferencesPanel(object):
    def setupUi(self, PreferencesPanel):
        PreferencesPanel.setObjectName("PreferencesPanel")
        PreferencesPanel.resize(443, 253)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(PreferencesPanel.sizePolicy().hasHeightForWidth())
        PreferencesPanel.setSizePolicy(sizePolicy)
        PreferencesPanel.setMinimumSize(QtCore.QSize(443, 253))
        PreferencesPanel.setMaximumSize(QtCore.QSize(443, 253))
        self.gridLayout = QtWidgets.QGridLayout(PreferencesPanel)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(PreferencesPanel)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.coreBox = QtWidgets.QSpinBox(PreferencesPanel)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.coreBox.sizePolicy().hasHeightForWidth())
        self.coreBox.setSizePolicy(sizePolicy)
        self.coreBox.setMinimum(1)
        self.coreBox.setMaximum(64)
        self.coreBox.setProperty("value", 4)
        self.coreBox.setObjectName("coreBox")
        self.gridLayout.addWidget(self.coreBox, 2, 1, 1, 1)
        self.applySettingsButton = QtWidgets.QPushButton(PreferencesPanel)
        self.applySettingsButton.setObjectName("applySettingsButton")
        self.gridLayout.addWidget(self.applySettingsButton, 4, 0, 1, 2)
        self.line_2 = QtWidgets.QFrame(PreferencesPanel)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout.addWidget(self.line_2, 3, 0, 1, 2)
        self.thresholdField = QtWidgets.QLineEdit(PreferencesPanel)
        self.thresholdField.setClearButtonEnabled(True)
        self.thresholdField.setObjectName("thresholdField")
        self.gridLayout.addWidget(self.thresholdField, 1, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(PreferencesPanel)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(PreferencesPanel)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 2, 0, 1, 1)
        self.nsPerFrameField = QtWidgets.QLineEdit(PreferencesPanel)
        self.nsPerFrameField.setClearButtonEnabled(True)
        self.nsPerFrameField.setObjectName("nsPerFrameField")
        self.gridLayout.addWidget(self.nsPerFrameField, 0, 1, 1, 1)

        self.retranslateUi(PreferencesPanel)
        QtCore.QMetaObject.connectSlotsByName(PreferencesPanel)
        PreferencesPanel.setTabOrder(self.nsPerFrameField, self.thresholdField)
        PreferencesPanel.setTabOrder(self.thresholdField, self.coreBox)
        PreferencesPanel.setTabOrder(self.coreBox, self.applySettingsButton)

    def retranslateUi(self, PreferencesPanel):
        _translate = QtCore.QCoreApplication.translate
        PreferencesPanel.setWindowTitle(_translate("PreferencesPanel", "Preferences"))
        self.label.setText(_translate("PreferencesPanel", "Nanoseconds/frame:"))
        self.applySettingsButton.setText(_translate("PreferencesPanel", "Apply"))
        self.thresholdField.setText(_translate("PreferencesPanel", "0.0"))
        self.label_2.setText(_translate("PreferencesPanel", "Score Threshold:"))
        self.label_4.setText(_translate("PreferencesPanel", "Cores:"))
        self.nsPerFrameField.setText(_translate("PreferencesPanel", "1.0"))

