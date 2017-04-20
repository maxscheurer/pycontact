from __future__ import print_function
import sip
import time
import os


from PyQt5.QtWidgets import QWidget, QApplication, QFileDialog
import MDAnalysis
import numpy as np

from .Plotters import ContactPlotter
from statistics_ui import *
from ..core.multi_accumulation import chunks
from ..core.Biochemistry import vdwRadius
from ..core.LogPool import *
from ..cy_modules import cy_gridsearch
from Dialogues import TopoTrajLoaderDialog
from ErrorBox import ErrorBox
from ErrorMessages import ErrorMessages
from ..core.Biochemistry import mean_score_of_contactArray, median_score_of_contactArray


class Statistics(QWidget, Ui_Statistics):
    def __init__(self, data, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)
        self.contacts = data

        self.labelNumFrames.setText(str(len(self.contacts[0].scoreArray)))
        self.labelTotalContacts.setText(str(len(self.contacts)))
        self.labelMeanScore.setText(str(mean_score_of_contactArray(self.contacts)))
        self.labelMedianScore.setText(str(median_score_of_contactArray(self.contacts)))

        self.contactPlotter = ContactPlotter(None, width=4, height=2, dpi=80)
        self.contactPlotter.plot_all_contacts_figure(self.contacts)
        # self.contactPlotter.plot_hbondNumber(self.contacts)
        self.plotGridLayout.addWidget(self.contactPlotter)
