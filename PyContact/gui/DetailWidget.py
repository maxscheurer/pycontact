from __future__ import print_function
import sip
import os

from PyQt5.QtWidgets import QWidget, QApplication, QFileDialog
import MDAnalysis
import numpy as np

from .Plotters import ContactPlotter
from detail_ui import *
from ..core.multi_accumulation import chunks
from ..core.Biochemistry import vdwRadius
from ..core.LogPool import *
from ..cy_modules import cy_gridsearch
from Dialogues import TopoTrajLoaderDialog
from ErrorBox import ErrorBox
from ErrorMessages import ErrorMessages
from ..core.Biochemistry import mean_score_of_contactArray, median_score_of_contactArray


class Detail(QWidget, Ui_Detail):
    def __init__(self, data, nsPerFrame, threshold, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)
        self.contact = data
        self.nsPerFrame = nsPerFrame
        self.threshold = threshold

        self.setWindowTitle(self.contact.title)
        self.labelTotalTime.setText(str(self.contact.total_time(self.nsPerFrame, self.threshold)))
        self.labelThreshold.setText(str(self.threshold))
        self.labelMedianScore.setText(str(self.contact.median_score()))
        self.labelMeanScore.setText(str(self.contact.mean_score()))
        self.labelMedianLifetime.setText(str(self.contact.median_life_time(self.nsPerFrame, self.threshold)))
        self.labelMeanLifetime.setText(str(self.contact.mean_life_time(self.nsPerFrame, self.threshold)))
        self.labelBackboneSidechainA.setText("%.2f/%.2f" % (self.contact.bb1, self.contact.sc1))
        self.labelBackboneSidechainB.setText("%.2f/%.2f" % (self.contact.bb2, self.contact.sc2))

        self.savePlotButton.clicked.connect(self.savePlot)
        self.plotButton.clicked.connect(self.plotAttribute)

        self.contactPlotter = ContactPlotter(None, width=4, height=2, dpi=70)
        self.contactPlotter.plot_contact_figure(self.contact, self.nsPerFrame)
        self.plotGridLayout.addWidget(self.contactPlotter)

    def plotAttribute(self):
        """Plots the selected attribute."""
        sip.delete(self.contactPlotter)
        self.contactPlotter = ContactPlotter(None, width=4, height=2, dpi=70)
        if self.attributeBox.currentText() == "Score":
            self.contactPlotter.plot_contact_figure(self.contact, self.nsPerFrame)
        self.plotGridLayout.addWidget(self.contactPlotter)

    def savePlot(self):
        """Saves the current plot."""
        fileName = QFileDialog.getSaveFileName(self, 'Save Path')
        path, file_extension = os.path.splitext(fileName[0])
        if file_extension == "":
            file_extension = "png"
        else:
            file_extension = file_extension[1:]
        try:
            self.contactPlotter.saveFigure(path, file_extension)
        except ValueError:
            box = ErrorBox("File format " + file_extension + " is not supported.\nPlease choose from eps, pdf, pgf,"
                                                             " png, ps, raw, rgba, svg, svgz. ")
            box.exec_()
