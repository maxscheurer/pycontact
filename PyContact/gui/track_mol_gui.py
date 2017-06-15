# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'track_mol_gui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_trackMoleculeView(object):
    def setupUi(self, trackMoleculeView):
        trackMoleculeView.setObjectName("trackMoleculeView")
        trackMoleculeView.resize(584, 443)
        self.gridLayout = QtWidgets.QGridLayout(trackMoleculeView)
        self.gridLayout.setObjectName("gridLayout")
        self.widget = QtWidgets.QWidget(trackMoleculeView)
        self.widget.setObjectName("widget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.trackTimelineView = QtWidgets.QWidget(self.widget)
        self.trackTimelineView.setMinimumSize(QtCore.QSize(0, 340))
        self.trackTimelineView.setObjectName("trackTimelineView")
        self.gridLayout_2.addWidget(self.trackTimelineView, 1, 0, 1, 1)
        self.trackSettingsView = QtWidgets.QWidget(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.trackSettingsView.sizePolicy().hasHeightForWidth())
        self.trackSettingsView.setSizePolicy(sizePolicy)
        self.trackSettingsView.setObjectName("trackSettingsView")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.trackSettingsView)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.sel1RadioButton = QtWidgets.QRadioButton(self.trackSettingsView)
        self.sel1RadioButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.sel1RadioButton.setChecked(True)
        self.sel1RadioButton.setObjectName("sel1RadioButton")
        self.horizontalLayout.addWidget(self.sel1RadioButton)
        self.sel2RadioButton = QtWidgets.QRadioButton(self.trackSettingsView)
        self.sel2RadioButton.setMaximumSize(QtCore.QSize(100, 16777215))
        self.sel2RadioButton.setObjectName("sel2RadioButton")
        self.horizontalLayout.addWidget(self.sel2RadioButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.label = QtWidgets.QLabel(self.trackSettingsView)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.mergeFrames = QtWidgets.QLineEdit(self.trackSettingsView)
        self.mergeFrames.setMaximumSize(QtCore.QSize(70, 16777215))
        self.mergeFrames.setObjectName("mergeFrames")
        self.horizontalLayout.addWidget(self.mergeFrames)
        self.runTrackingButton = QtWidgets.QPushButton(self.trackSettingsView)
        self.runTrackingButton.setObjectName("runTrackingButton")
        self.horizontalLayout.addWidget(self.runTrackingButton)
        self.gridLayout_2.addWidget(self.trackSettingsView, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 1)

        self.retranslateUi(trackMoleculeView)
        QtCore.QMetaObject.connectSlotsByName(trackMoleculeView)

    def retranslateUi(self, trackMoleculeView):
        _translate = QtCore.QCoreApplication.translate
        trackMoleculeView.setWindowTitle(_translate("trackMoleculeView", "Form"))
        self.sel1RadioButton.setText(_translate("trackMoleculeView", "Sel. 1"))
        self.sel2RadioButton.setText(_translate("trackMoleculeView", "Sel. 2"))
        self.label.setText(_translate("trackMoleculeView", "Merge frames:"))
        self.mergeFrames.setText(_translate("trackMoleculeView", "1"))
        self.runTrackingButton.setText(_translate("trackMoleculeView", "Run Tracking"))

