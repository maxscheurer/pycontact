import sip
import os

from PyQt5.QtWidgets import QWidget, QFileDialog
from PyQt5.QtGui import QIntValidator

from .Plotters import ContactPlotter
from .statistics_ui import *
from ..core.LogPool import *
from .ErrorBox import ErrorBox
from ..core.Biochemistry import mean_score_of_contactArray, median_score_of_contactArray


class Statistics(QWidget, Ui_Statistics):
    def __init__(self, data, nspf, parent=None):
        super(QWidget, self).__init__(parent)
        self.setupUi(self)
        self.contacts = data
        self.nsPerFrame = nspf

        self.labelNumFrames.setText(str(len(self.contacts[0].scoreArray)))
        self.labelTotalContacts.setText(str(len(self.contacts)))
        self.labelMeanScore.setText(str(mean_score_of_contactArray(self.contacts)))
        self.labelMedianScore.setText(str(median_score_of_contactArray(self.contacts)))
        self.savePlotButton.clicked.connect(self.savePlot)
        self.plotButton.clicked.connect(self.plotAttribute)

        posIntValidator = QIntValidator()
        posIntValidator.setBottom(3)
        self.smoothStrideField.setValidator(posIntValidator)

        self.contactPlotter = ContactPlotter(None, width=4, height=2, dpi=70)
        self.contactPlotter.plot_all_contacts_figure(self.contacts, 0, self.nsPerFrame)
        self.plotGridLayout.addWidget(self.contactPlotter)

    def plotAttribute(self):
        """Plots the selected attribute."""
        sip.delete(self.contactPlotter)
        self.contactPlotter = ContactPlotter(None, width=4, height=2, dpi=70)
        smoothOn = self.smoothCheckbox.isChecked()
        smooth = 0
        limit = 5
        if smoothOn:
            smooth = int(self.smoothStrideField.text())
            if smooth < limit:
                smooth = limit
            if smooth % 2 == 0:
                smooth += 1
            self.smoothStrideField.setText(str(smooth))
        if self.attributeBox.currentText() == "Score":
            self.contactPlotter.plot_all_contacts_figure(self.contacts, smooth, self.nsPerFrame)
        elif self.attributeBox.currentText() == "hbond number":
            self.contactPlotter.plot_hbondNumber(self.contacts, smooth, self.nsPerFrame)
        self.plotGridLayout.addWidget(self.contactPlotter)

    def savePlot(self):
        """Saves the current plot."""
        fileName = QFileDialog.getSaveFileName(self, 'Save Path')
        if fileName == "":
            return
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
