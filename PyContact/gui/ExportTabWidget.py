import sip
import os
from PyQt5.QtWidgets import QTabWidget, QWidget, QGridLayout, QLabel, QPushButton, QComboBox, QLineEdit, QCheckBox, \
    QFileDialog
from PyQt5.QtCore import pyqtSignal

from Plotters import *
from ErrorBox import ErrorBox
from ErrorMessages import ErrorMessages


class ExportTabWidget(QTabWidget):
    """Widget for data exporting purposes."""

    valueUpdated = pyqtSignal(str, str)

    def __init__(self, parent=None):
        super(ExportTabWidget, self).__init__(parent)
        self.checkboxes = []
        self.keys = []
        self.checkboxdict = []
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()
        self.tab5 = QWidget()
        self.grid1 = QGridLayout()
        self.grid2 = QGridLayout()
        self.grid3 = QGridLayout()
        self.grid4 = QGridLayout()

        self.addTab(self.tab1, "View")
        self.addTab(self.tab2, "Histogram")
        self.addTab(self.tab3, "Contact Map")
        self.addTab(self.tab4, "VMD")
        self.addTab(self.tab5, "Plain Text")
        self.tab1UI()
        self.tab2UI()
        self.tab3UI()
        self.tab4UI()
        self.tab5UI()
        self.setWindowTitle("Export")
        self.contacts = []
        self.threshold = 0
        self.nsPerFrame = 0
        self.map1 = None
        self.map2 = None
        self.label1 = ""
        self.label2 = ""
        self.psf = ""
        self.dcd = ""

    def setThresholdAndNsPerFrame(self, currentThreshold, currentNsPerFrame):
        self.threshold = currentThreshold
        self.nsPerFrame = currentNsPerFrame

    def tab1UI(self):
        """Tab where the current view can be exported to file."""
        grid = QGridLayout()
        self.tab1.setLayout(grid)

        self.tab1.exportLabel = QLabel("Export current view: ")

        self.tab1.saveButton = QPushButton("Export")
        self.tab1.saveButton.setAutoDefault(False)
        self.tab1.saveButton.clicked.connect(self.pushSave)

        self.tab1.formatBox = QComboBox()
        self.tab1.formatBox.addItem("PNG")
        self.tab1.formatBox.addItem("SVG")

        grid.addWidget(self.tab1.exportLabel, 0, 0)
        grid.addWidget(self.tab1.saveButton, 0, 1)

        grid.addWidget(self.tab1.formatBox, 2, 0)

    def tab2UI(self):
        """Tab where analyzed data can be visualized and exported as histograms."""
        self.tab2.setLayout(self.grid1)

        self.tab2.histPlot = HistPlotter(None, width=8, height=5, dpi=60)
        self.grid1.addWidget(self.tab2.histPlot, 3, 0, 1, 4)

        self.tab2.histTypeBox = QComboBox()
        self.tab2.histTypeBox.addItem("General Histogram")
        self.tab2.histTypeBox.addItem("Bin per Contact")

        self.grid1.addWidget(self.tab2.histTypeBox, 0, 3)

        self.tab2.attributeBox = QComboBox()
        self.tab2.attributeBox.addItem("Mean Score")
        self.tab2.attributeBox.addItem("Median Score")
        self.tab2.attributeBox.addItem("Mean Lifetime")
        self.tab2.attributeBox.addItem("Median Lifetime")
        self.tab2.attributeBox.addItem("Hbond percentage")

        self.grid1.addWidget(self.tab2.attributeBox, 1, 3)

        self.tab2.plotButton = QPushButton("Show Preview")
        self.tab2.plotButton.setAutoDefault(False)
        self.tab2.plotButton.clicked.connect(self.pushPlot)
        self.grid1.addWidget(self.tab2.plotButton, 0, 0, 1, 3)

        self.tab2.saveButton = QPushButton("Save Histogram")
        self.tab2.saveButton.setAutoDefault(False)
        self.tab2.saveButton.clicked.connect(self.saveHist)
        self.grid1.addWidget(self.tab2.saveButton, 1, 0, 1, 1)

        self.tab2.formatLabel = QLabel("Format: ")
        self.grid1.addWidget(self.tab2.formatLabel, 1, 1)

        self.tab2.xTicksFontSizeLabel = QLabel("bin per contact font size: ")
        self.grid1.addWidget(self.tab2.xTicksFontSizeLabel, 2, 0)

        self.tab2.xTicksFontSizeField = QLineEdit("11")
        self.grid1.addWidget(self.tab2.xTicksFontSizeField, 2, 1)

        self.tab2.formatBox = QComboBox()
        self.tab2.formatBox.addItem("pdf")
        self.tab2.formatBox.addItem("png")
        self.tab2.formatBox.addItem("svg")
        self.tab2.formatBox.addItem("eps")
        self.grid1.addWidget(self.tab2.formatBox, 1, 2)

    def tab3UI(self):
        """Tab where the contact map can be generated and exported."""
        self.tab3.setLayout(self.grid2)

        self.tab3.mapPlot = AnimateMapPlotter(None, width=8, height=5, dpi=60)
        self.grid2.addWidget(self.tab3.mapPlot, 3, 0, 1, 4)

        self.tab3.plotMapButton = QPushButton("Show Preview")
        self.tab3.plotMapButton.setAutoDefault(False)
        self.tab3.plotMapButton.clicked.connect(self.pushMapPlot)
        self.grid2.addWidget(self.tab3.plotMapButton, 0, 0, 1, 1)

        self.tab3.formatBox = QComboBox()
        self.tab3.formatBox.addItem("pdf")
        self.tab3.formatBox.addItem("png")
        self.tab3.formatBox.addItem("svg")
        self.tab3.formatBox.addItem("eps")
        self.grid2.addWidget(self.tab3.formatBox, 0, 2, 1, 1)

        self.tab3.saveButton = QPushButton("Save Map")
        self.tab3.saveButton.setAutoDefault(False)
        self.tab3.saveButton.clicked.connect(self.saveMap)
        self.grid2.addWidget(self.tab3.saveButton, 0, 3, 1, 1)

        self.tab3.attributeBox = QComboBox()
        self.tab3.attributeBox.addItem("Mean Score")
        self.tab3.attributeBox.addItem("Median Score")
        self.tab3.attributeBox.addItem("Mean Lifetime")
        self.tab3.attributeBox.addItem("Median Lifetime")
        self.tab3.attributeBox.addItem("Hbond percentage")

        self.grid2.addWidget(self.tab3.attributeBox, 0, 1)

    def tab4UI(self):
        """Tab where selected data can be visualized in VMD via a Tcl script generation."""
        self.tab4.setLayout(self.grid3)

        label = QLabel("Split selections for each contact")
        self.tab4.splitVisCheckbox = QCheckBox()
        self.grid3.addWidget(label, 0, 0)
        self.grid3.addWidget(self.tab4.splitVisCheckbox, 0, 1)

        self.tab4.additionalText1 = QLineEdit()
        self.tab4.additionalText2 = QLineEdit()
        additionalTextLabel1 = QLabel("Additional seltext for selection 1")
        additionalTextLabel2 = QLabel("Additional seltext for selection 2")

        self.tab4.button = QPushButton("Create tcl script")
        self.tab4.button.clicked.connect(self.createTclScriptVis)
        self.grid3.addWidget(self.tab4.button, 3, 0, 1, 2)

        self.grid3.addWidget(additionalTextLabel1, 1, 0)
        self.grid3.addWidget(additionalTextLabel2, 2, 0)
        self.grid3.addWidget(self.tab4.additionalText1, 1, 1)
        self.grid3.addWidget(self.tab4.additionalText2, 2, 1)

    def tab5UI(self):
        """Tab where the raw data can be exported to a text file."""
        self.tab5.setLayout(self.grid4)

        self.checkboxdict = {"mean_score": "Mean Score",
                             "hbond_percentage": "HBond Percentage", "median_score": "Median Score",
                             "contactTypeAsShortcut": "Contact Type",
                             "getScoreArray": "Score List", "hbondFramesScan": "Hydrogen Bond Frames"}
        # let's define the key order ourselves
        self.keys = ["contactTypeAsShortcut", "mean_score", "median_score", "hbond_percentage", "getScoreArray", "hbondFramesScan"]

        propertyLabel = QLabel("Select properties to export")
        self.grid4.addWidget(propertyLabel, 0, 0)

        startLine = 1
        for box in self.keys:
            checkBox = QCheckBox()
            checkBox.setChecked(True)
            boxLabel = QLabel(self.checkboxdict[box])
            self.grid4.addWidget(boxLabel, startLine, 0)
            self.grid4.addWidget(checkBox, startLine, 1)
            self.checkboxes.append(checkBox)
            startLine += 1

        self.tab5.button = QPushButton("Export to text")
        self.tab5.button.clicked.connect(self.saveText)
        self.grid4.addWidget(self.tab5.button, startLine, 5)

    def saveText(self):
        """Executes the conversion of the analysis results to raw text data."""
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            boxIndex = 0
            requestedParameters = []
            for box in self.checkboxes:
                if box.isChecked():
                    requestedParameters.append(self.keys[boxIndex])
                boxIndex += 1

            tableHeadings = []
            for par in requestedParameters:
                tableHeadings.append(self.checkboxdict[par])

            f = open(path + ".txt", "w")
            row_format = " {:>20} " * (len(requestedParameters) + 1)
            f.write(row_format.format("", *tableHeadings))
            f.write("\n")

            for c in self.contacts:
                currentContactProperties = []
                for p in requestedParameters:
                    exec("propertyToAdd = c." + p + "()")
                    if isinstance(propertyToAdd, float):
                        propertyToAdd = "{0:.3f}".format(propertyToAdd)
                    currentContactProperties.append(propertyToAdd)
                f.write(row_format.format(c.human_readable_title(), *currentContactProperties))
                f.write("\n")
            f.close()

    def saveHist(self):
        """Saves the histogram to the picked file path."""
        self.plotHist()

        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            if file_extension == "":
                file_extension = "." + self.tab2.formatBox.currentText().lower()
            path += file_extension
            self.tab2.histPlot.saveFigure(path, self.tab2.formatBox.currentText())

    def saveMap(self):
        """Saves the contact map to the picked file path."""
        self.plotMap()

        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            self.tab3.mapPlot.saveFigure(path, self.tab3.formatBox.currentText())

    def pushPlot(self):
        """Triggers the histogram plotter."""
        self.plotHist()

    def pushMapPlot(self):
        """Triggers the contact map plotter."""
        self.plotMap()

    def plotHist(self):
        """Plots the histogram."""
        sip.delete(self.tab2.histPlot)
        self.tab2.histPlot = HistPlotter(None, width=8, height=5, dpi=60)
        self.grid1.addWidget(self.tab2.histPlot, 3, 0, 1, 4)
        if self.tab2.histTypeBox.currentText() == "General Histogram":
            self.tab2.histPlot.plotGeneralHist(self.contacts, self.tab2.attributeBox.currentText(), self.threshold,
                                               self.nsPerFrame)
        elif self.tab2.histTypeBox.currentText() == "Bin per Contact":
            self.tab2.histPlot.plotContactHist(self.contacts, self.tab2.attributeBox.currentText(), self.threshold,
                                               self.nsPerFrame, int(self.tab2.xTicksFontSizeField.text()))

        self.tab2.histPlot.update()

    def plotMap(self):
        """Plots the contact map."""
        sip.delete(self.tab3.mapPlot)
        self.tab3.mapPlot = MapPlotter(None, width=8, height=5, dpi=60)
        self.grid2.addWidget(self.tab3.mapPlot, 3, 0, 1, 4)
        if self.map1 is None or self.map2 is None or self.contacts is None or len(self.contacts) == 0:
            box = ErrorBox(ErrorMessages.RESID_REQUIRED)
            box.exec_()
            return
        res = self.tab3.mapPlot.plotMap(self.contacts, self.map1, self.map2, self.label1, self.label2,
                                        self.tab3.attributeBox.currentText(), self.threshold, self.nsPerFrame)
        if res == -1:
            box = ErrorBox(ErrorMessages.RESID_REQUIRED)
            box.exec_()
        self.tab3.mapPlot.update()

    def pushSave(self):
        """Saves the current view."""
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        path, file_extension = os.path.splitext(fileName[0])
        if file_extension == "":
            file_extension = "." + self.tab1.formatBox.currentText().lower()
        path += file_extension
        self.valueUpdated.emit(path, self.tab1.formatBox.currentText())

    def setContacts(self, currentContacts):
        self.contacts = currentContacts

    def setFilePaths(self, *argv):
        """Sets the current trajectory paths from the main view."""
        self.psf = argv[0][0]
        self.dcd = argv[0][1]

    def setMaps(self, map1, map2):
        self.map1 = map1
        self.map2 = map2

# labels for the contact map
    def setMapLabels(self, label1, label2):
        self.label1 = label1
        self.label2 = label2

    def createTclScriptVis(self):
        """Creates the Tcl script for VMD visualization of the selections."""
        if len(self.contacts) == 0:
            box = ErrorBox(ErrorMessages.NOCONTACTS)
            box.exec_()
            return

        fileName = QFileDialog.getSaveFileName(self, 'Save Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            if file_extension != ".tcl":
                file_extension = ".tcl"

            path += file_extension

            f = open(path, 'w')
            f.write('mol new %s \n' % self.psf)
            f.write('mol addfile %s \n' % self.dcd)
            f.write('mol delrep 0 top \n')
            f.write('mol representation NewCartoon \n')
            f.write('mol Color ColorID 3 \n')
            f.write('mol selection {all} \n')
            f.write('mol addrep top \n')

            if self.tab4.splitVisCheckbox.isChecked():
                for cont in self.contacts:
                    currentSel1 = []
                    index = 0
                    for item in cont.key1:
                        if item != "none":
                            currentSel1.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                        index += 1
                    currentSel1String = " and ".join(currentSel1)
                    currentSel2 = []
                    index = 0
                    for item in cont.key2:
                        if item != "none":
                            currentSel2.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                        index += 1
                    currentSel2String = " and ".join(currentSel2)
                    add1 = ("" if self.tab4.additionalText1.text() == "" else (" and " + self.tab4.additionalText1.text()))
                    add2 = ("" if self.tab4.additionalText2.text() == "" else (" and " + self.tab4.additionalText2.text()))
                    sel = "("+currentSel1String + add1 + ") or (" + currentSel2String + add2 + ")"
                    f.write('mol representation Licorice \n')
                    f.write('mol Color Name \n')
                    f.write('mol selection {%s} \n' % sel)
                    f.write('mol addrep top \n')
            else:
                total = []
                for cont in self.contacts:
                    currentSel1 = []
                    index = 0
                    for item in cont.key1:
                        if item != "none":
                            currentSel1.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                        index += 1
                    currentSel1String = " and ".join(currentSel1)
                    currentSel2 = []
                    index = 0
                    for item in cont.key2:
                        if item != "none":
                            currentSel2.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                        index += 1
                    currentSel2String = " and ".join(currentSel2)
                    add1 = ("" if self.tab4.additionalText1.text() == "" else (" and " + self.tab4.additionalText1.text()))
                    add2 = ("" if self.tab4.additionalText2.text() == "" else (" and " + self.tab4.additionalText2.text()))
                    sel = "(" + currentSel1String + add1 + ") or (" + currentSel2String + add2 + ")"
                    total.append(sel)
                seltext = " or ".join(total)
                f.write('mol representation Licorice \n')
                f.write('mol Color Name \n')
                f.write('mol selection {%s} \n' % seltext)
                f.write('mol addrep top \n')
            f.close()
