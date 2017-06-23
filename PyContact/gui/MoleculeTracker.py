from __future__ import print_function
import sip, os

from PyQt5.QtWidgets import QWidget, QRadioButton, QApplication, QFileDialog, QLabel
from PyQt5.QtGui import (QColor, QPainter, QFont, QPixmap, QPaintEvent)
from PyQt5.QtCore import QSize, QRect
from PyQt5.QtCore import pyqtSignal, QObject
import numpy as np

from track_mol_gui import *
from ErrorMessages import ErrorMessages
from ErrorBox import ErrorBox
from TrackCanvas import TrackCanvas


class MoleculeTracker(QWidget, Ui_trackMoleculeView):
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("MoleculeTracker")
        self.contactAnalyzer = None
        self.runTrackingButton.clicked.connect(self.runTracking)
        self.textExportButton.clicked.connect(self.exportTextFile)
        self.labelPainter = TrackCanvas()
        self.scrollArea.setWidget(self.labelPainter)
        self.trackingResult = None

    def setContactAnalyzer(self, conAnalyzer):
        self.contactAnalyzer = conAnalyzer

    def clean(self):
        self.contactAnalyzer = None

    def exportTextFile(self):
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            if file_extension == "":
                file_extension = ".dat"

            f = open(path + file_extension, "w")
            for i in range(len(self.trackingResult)):
                cont = self.trackingResult[i]
                line2print = [(x[0]+" ("+str(x[1])+") ") for x in cont[:int(self.maxContactsPerFrame.text())]]
                f.write(str(i) + "\t" + " ".join(line2print) + "\n")
            f.close()

    def runTracking(self):
        if self.contactAnalyzer is None:
            box = ErrorBox(ErrorMessages.NODATA_PROMPTLOAD)
            box.exec_()
            return
        print("Executing molecule tracking.")
        selectionIndex = 0
        if self.sel1RadioButton.isChecked():
            selectionIndex = 1
        elif self.sel2RadioButton.isChecked():
            selectionIndex = 2
        # print("Selection: ", selectionIndex)
        frameMerge = int(self.mergeFrames.text())
        # print("Merge: ", frameMerge)
        self.trackingResult = self.contactAnalyzer.runMoleculeTracking(selectionIndex, [0, 0, 1, 1, 0])
        # print("Analysis finished", result)
        # self.labelPainter.contacts = result
        # self.labelPainter.maximalContactsPerRow = int(self.maxContactsPerFrameField.text())
        # self.labelPainter.draw_labels()
        # self.labelPainter.repaint()
        # self.labelPainter.update()
