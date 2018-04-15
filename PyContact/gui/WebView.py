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

    def setUrl(self):
        self.webEngineView.setUrl(QUrl("https://youtube.com/"))

    def setData(self, data):
        self.webEngineView.setHtml(data, QUrl("./"))

    def plotBokeh(self):
        n = 500
        x = 2 + 2*np.random.standard_normal(n)
        y = 2 + 2*np.random.standard_normal(n)

        p = figure(title="Hexbin for 500 points", match_aspect=True,
                   tools="wheel_zoom,reset", background_fill_color='#440154')
        p.grid.visible = False

        r, bins = p.hexbin(x, y, size=0.5, hover_color="pink", hover_alpha=0.8)

        p.circle(x, y, color="white", size=1)

        hover = HoverTool(tooltips=[("count", "@c"), ("(q,r)", "(@q, @r)")],
                          mode="mouse", point_policy="follow_mouse", renderers=[r])

        p.add_tools(hover)

        html = file_html(p, CDN, "Test Plot")
        self.setData(html)
