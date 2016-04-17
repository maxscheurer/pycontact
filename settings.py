# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'settings.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_settingsWidget(object):
    def setupUi(self, settingsWidget):
        settingsWidget.setObjectName("settingsWidget")
        settingsWidget.resize(574, 447)
        self.filterTab = QtWidgets.QWidget()
        self.filterTab.setObjectName("filterTab")
        settingsWidget.addTab(self.filterTab, "")
        self.settingsTab = QtWidgets.QWidget()
        self.settingsTab.setObjectName("settingsTab")
        settingsWidget.addTab(self.settingsTab, "")

        self.retranslateUi(settingsWidget)
        settingsWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(settingsWidget)

    def retranslateUi(self, settingsWidget):
        _translate = QtCore.QCoreApplication.translate
        settingsWidget.setWindowTitle(_translate("settingsWidget", "Filters and Settings"))
        settingsWidget.setTabText(settingsWidget.indexOf(self.filterTab), _translate("settingsWidget", "Filters"))
        settingsWidget.setTabText(settingsWidget.indexOf(self.settingsTab), _translate("settingsWidget", "Settings"))

