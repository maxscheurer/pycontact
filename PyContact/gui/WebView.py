from __future__ import print_function

from PyQt5.QtWidgets import QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
import numpy as np
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool

from .webview_gui import *



class WebView(QWidget, Ui_WebView):
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)
        self.gasButton.clicked.connect(self.setUrl)
        self.trajectory = None

    def setUrl(self):
        self.webEngineView.setUrl(QUrl("https://youtube.com/"))

    def setData(self, data):
        self.webEngineView.setHtml(data, QUrl("./"))

    def setAccumulatedTrajectory(self, trajectory):
        self.trajectory = trajectory

    def plotBokeh(self):
        
