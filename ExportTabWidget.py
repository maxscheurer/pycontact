from PyQt5.QtWidgets import QTabWidget, QWidget, QGridLayout, QLabel, QPushButton, QComboBox, QLineEdit
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from Plotters import *
class ExportTabWidget(QTabWidget):

    valueUpdated = pyqtSignal(str, str)

    def __init__(self, parent = None):
       super(ExportTabWidget, self).__init__(parent)

       self.tab1 = QWidget()
       self.tab2 = QWidget()

       self.addTab(self.tab1, "View")
       self.addTab(self.tab2, "Histogram")
       self.tab1UI()
       self.tab2UI()
       self.setWindowTitle("Export")
       self.contacts = []
       self.threshold = 0
       self.nsPerFrame = 0

    def setThresholdAndNsPerFrame(self, currentThreshold, currentNsPerFrame):
        self.threshold = currentThreshold
        self.nsPerFrame = currentNsPerFrame

    def tab1UI(self):
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
        self.grid1 = QGridLayout()
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

    def saveHist(self):
        self.plotHist()

        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            self.tab2.histPlot.saveFigure(path, self.tab2.formatBox.currentText())

    def pushPlot(self):
        self.plotHist()

    def plotHist(self):
        sip.delete(self.tab2.histPlot)
        self.tab2.histPlot = HistPlotter(None, width=8, height=5, dpi=60)
        self.grid1.addWidget(self.tab2.histPlot, 3, 0, 1, 4)
        if self.tab2.histTypeBox.currentText() == "General Histogram":
            self.tab2.histPlot.plotGeneralHist(self.contacts, self.tab2.attributeBox.currentText(), self.threshold, self.nsPerFrame)
        elif self.tab2.histTypeBox.currentText() == "Bin per Contact":
            self.tab2.histPlot.plotContactHist(self.contacts, self.tab2.attributeBox.currentText(), self.threshold, self.nsPerFrame,int(self.tab2.xTicksFontSizeField.text()))

        self.tab2.histPlot.update()

    def pushSave(self):
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        self.valueUpdated.emit(fileName[0], self.tab1.formatBox.currentText())

    def setContacts(self, currentContacts):
        self.contacts = currentContacts
