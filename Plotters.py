import sys, sip, copy
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget, QDialog, QTabWidget, QButtonGroup,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QVBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton, QRadioButton,
                             QFrame, QApplication, QSizePolicy)
from PyQt5.QtSvg import QSvgGenerator
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore"); 
    import matplotlib.pyplot as plt
import gui
from settings import *
from biochemistry import *
from inputreader import *
from filters import *
from functools import partial

class MplPlotter(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class ContactPlotter(MplPlotter):
    """Simple canvas with a sine plot."""

    def plot_contact_figure(self, contact):
        self.axes.plot(contact.scoreArray)
        self.axes.set_xlabel("frame")
        self.axes.set_ylabel("score")

    def plot_all_contacts_figure(self, contacts):
        values = []

        for frame in range(len(contacts[0].scoreArray)):
            current = 0
            for c in contacts:
                current += c.scoreArray[frame]
            values.append(current)
        self.axes.plot(values)
        self.axes.set_xlabel("frame")
        self.axes.set_ylabel("score")


class HistPlotter(MplPlotter):
    """Simple canvas with an histogram plot."""

    def plotGeneralHist(self, currentContacts, attribute, threshold, nsPerFrame):
        values = []

        if attribute == "Mean Score":
            for c in currentContacts:
                values.append(c.mean_score())
        elif attribute == "Median Score":
            for c in currentContacts:
                values.append(c.median_score())
        elif attribute == "Mean Lifetime":
            for c in currentContacts:
                values.append(c.mean_life_time(nsPerFrame, threshold))
        elif attribute == "Median Lifetime":
            for c in currentContacts:
                values.append(c.median_life_time(nsPerFrame, threshold))

        valuesNp = np.array(values, dtype = float)
        self.axes.hist(valuesNp, bins=20)
        self.fig.subplots_adjust(bottom=0.2, top=0.95, left=0.15, right=0.85)

    def plotContactHist(self, currentContacts, attribute, threshold, nsPerFrame):
        values = []
        titles = []

        if attribute == "Mean Score":
            for c in currentContacts:
                values.append(c.mean_score())
                titles.append(c.title)
        elif attribute == "Median Score":
            for c in currentContacts:
                values.append(c.median_score())
                titles.append(c.title)
        elif attribute == "Mean Lifetime":
            for c in currentContacts:
                values.append(c.mean_life_time(nsPerFrame, threshold))
                titles.append(c.title)
        elif attribute == "Median Lifetime":
            for c in currentContacts:
                values.append(c.median_life_time(nsPerFrame, threshold))
                titles.append(c.title)

        valuesNp = np.array(values, dtype = float)
        titlesNp = np.array(titles, dtype = str)

        x = range(len(currentContacts))
        h = self.axes.bar(x, valuesNp, color="red")
        xticks_pos = [0.7071 * patch.get_width() + patch.get_xy()[0] for patch in h]
        self.axes.set_xticklabels(titlesNp, ha='right', size=8, rotation=45)
        self.axes.set_xticks(xticks_pos)
        self.fig.subplots_adjust(bottom=0.2, top=0.95, left=0.1, right=0.9)

    def saveFigure(self, path, outputFormat):
        self.fig.savefig(path + "." + outputFormat, format=outputFormat)


class SimplePlotter(MplPlotter):
    """Simple canvas with a sine plot."""
    def plot(self, x, y):
        self.axes.plot(x, y)
        self.axes.set_xlabel("x")
        self.axes.set_ylabel("f(x)")
        self.axes.xaxis.set_label_position('top')