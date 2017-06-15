from __future__ import print_function
import sip

from PyQt5.QtWidgets import QWidget, QRadioButton, QApplication, QFileDialog

from track_mol_gui import *
from ErrorMessages import ErrorMessages
from ErrorBox import ErrorBox


class MoleculeTracker(QWidget, Ui_trackMoleculeView):
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("MoleculeTracker")
        self.contactAnalyzer = None
        self.runTrackingButton.clicked.connect(self.runTracking)

    def setContactAnalyzer(self, conAnalyzer):
        self.contactAnalyzer = conAnalyzer

    def clean(self):
        self.contactAnalyzer = None

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
        print("Selection: ", selectionIndex)
        frameMerge = int(self.mergeFrames.text())
        print("Merge: ", frameMerge)
        result = self.contactAnalyzer.runMoleculeTracking(selectionIndex, [0, 0, 1, 1, 0])
