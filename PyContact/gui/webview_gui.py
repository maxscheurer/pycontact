# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'webview.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_WebView(object):
    def setupUi(self, WebView):
        WebView.setObjectName("WebView")
        WebView.resize(1200, 700)
        self.gridLayout = QtWidgets.QGridLayout(WebView)
        self.gridLayout.setObjectName("gridLayout")
        self.webEngineView = QtWebEngineWidgets.QWebEngineView(WebView)
        self.webEngineView.setUrl(QtCore.QUrl("https://www.google.de/"))
        self.webEngineView.setObjectName("webEngineView")
        self.gridLayout.addWidget(self.webEngineView, 0, 0, 1, 1)
        self.gasButton = QtWidgets.QPushButton(WebView)
        self.gasButton.setObjectName("gasButton")
        self.gridLayout.addWidget(self.gasButton, 1, 0, 1, 1)

        self.retranslateUi(WebView)
        QtCore.QMetaObject.connectSlotsByName(WebView)

    def retranslateUi(self, WebView):
        _translate = QtCore.QCoreApplication.translate
        WebView.setWindowTitle(_translate("WebView", "Surface Areas"))
        self.gasButton.setText(_translate("WebView", "Gas!"))

from PyQt5 import QtWebEngineWidgets
