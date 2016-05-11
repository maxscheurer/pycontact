'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
from Canvas import *
from Plotters import *
from mdanalysis import *
from PyQt5.QtWidgets import *
import matplotlib.pyplot as plt

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.mlab import bivariate_normal
from mpl_toolkits.mplot3d import Axes3D

class MainWindow(QMainWindow, gui.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.contacts = []
        self.setupUi(self)

        self.setWindowTitle("pyContact")
        self.mergeSlider.valueChanged.connect(self.mergeValueChanged)

        # painter contains both labels and frame boxes for drawing
        self.painter = Canvas()
        self.scrollArea.setWidget(self.painter)
        self.actionOpen.triggered.connect(self.pushOpen)
        self.actionExport.triggered.connect(self.pushExport)
        # settings and filters
        self.settingsView = SettingsTabWidget()
        self.settingsView.applySettingsButton.clicked.connect(self.updateSettings)
        self.settingsView.applyFilterButton.clicked.connect(self.updateFilters)

        #alpha slider for color
        self.alphaSlider.setValue(50)
        self.alphaSlider.valueChanged.connect(self.alphaValueChanged)

        self.statisticsButton.clicked.connect(self.showStatistics)

        self.openPreferencesButton.clicked.connect(self.openPrefs)

        #color picker
        self.settingsView.pickColorButton.clicked.connect(self.showColorPicker)
        self.customColor = QColor(230, 50, 0)
        self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" % self.customColor.name())

        #visualize button
        self.visButton.clicked.connect(self.visualize)

        #button group for weight functions
        self.functionButtonGroup = QButtonGroup()
        self.currentFunctionType = FunctionType.sigmoid
        self.functionButtonGroup.buttonClicked[int].connect(self.showFunctionSettings)
        self.functionButtonGroup.addButton(self.settingsView.sigmoidRadioButton, 0)
        self.functionButtonGroup.addButton(self.settingsView.rectRadioButton, 1)
        self.functionButtonGroup.addButton(self.settingsView.linRadioButton, 2)
        self.setupFunctionBox()
        self.showFunctionSettings(FunctionType.sigmoid)

        #apply color button
        self.settingsView.applyColorButton.clicked.connect(self.updateColors)
        self.colorScheme = ColorScheme.bbsc

        map1 = [0, 0, 0, 1, 1, 0]
        map2 = [0, 0, 0, 1, 1, 0]
        self.map1 = map1
        self.map2 = map2
        self.psf = "rpn11_ubq_interface-ionized.psf"
        self.dcd = "short.dcd"
        contactResults = analyze_psf_dcd(self.psf, self.dcd, 5.0, 2.5, 120, "segid RN11","segid UBQ")
        self.contacts= analyze_contactResultsWithMaps(contactResults, map1, map2)

        # map1 = [0, 0, 0, 1, 1, 0]
        # map2 = [0, 0, 0, 1, 1, 0]
        # contactResults = analyze_psf_dcd("membraneAligned-helixZ90X0-ionized.psf", "arf_hmmm_100ns_100f_comet.dcd", 5.0, 2.5, 120, "segid PROT","segid MEMB")
        # self.contacts= analyze_contactResultsWithMaps(contactResults, map1, map2)

        # testing shit
        maxresids1 = []
        maxresids2 = []
        for cont in self.contacts:
            cont.determineBackboneSidechainType()
            maxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
            maxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))

        # Generate some test data
        x = np.arange(1,np.max(maxresids1)+2)
        y = np.arange(1,np.max(maxresids2)+2)
        data = np.zeros((len(x), len(y)))
        for cont in self.contacts:
            r1 = int(cont.key1[AccumulationMapIndex.resid])
            r2 = int(cont.key2[AccumulationMapIndex.resid])
            hbonds = cont.hbondFramesScan()
            count = np.count_nonzero(hbonds)
            if count > 0:
                print cont.title + " contains " + str(count) + " hbonds in total"
            # data[r1,r2] = cont.total_time(1, 0)
            data[r1, r2] = cont.mean_score()
        fig = plt.figure()
        ax = fig.add_subplot(111)

        cax = ax.matshow(data, cmap=cm.coolwarm, aspect='equal')
        fig.tight_layout()
        fig.colorbar(cax)
        # plt.show()
        self.updateSettings()
        self.updateFilters()

    def setupFunctionBox(self):
        # sig
        self.sigX0Label = QLabel("x0: ", self.settingsView)
        self.sigX0Field = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.sigX0Label, 1, 0)
        self.settingsView.functionGridLayout.addWidget(self.sigX0Field, 1, 1)

        self.sigLLabel = QLabel("L: ", self.settingsView)
        self.sigLField = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.sigLLabel, 1, 2)
        self.settingsView.functionGridLayout.addWidget(self.sigLField, 1, 3)

        self.sigKLabel = QLabel("k: ", self.settingsView)
        self.sigKField = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.sigKLabel, 2, 0)
        self.settingsView.functionGridLayout.addWidget(self.sigKField, 2, 1)

        self.sigY0Label = QLabel("y0: ", self.settingsView)
        self.sigY0Field = QLineEdit("0", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.sigY0Label, 2, 2)
        self.settingsView.functionGridLayout.addWidget(self.sigY0Field, 2, 3)

        # rect
        self.rectX0Label = QLabel("x0: ", self.settingsView)
        self.rectX0Field = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.rectX0Label, 1, 0)
        self.settingsView.functionGridLayout.addWidget(self.rectX0Field, 1, 1)

        self.rectX1Label = QLabel("x1: ", self.settingsView)
        self.rectX1Field = QLineEdit("2", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.rectX1Label, 1, 2)
        self.settingsView.functionGridLayout.addWidget(self.rectX1Field, 1, 3)

        self.rectHLabel = QLabel("h: ", self.settingsView)
        self.rectHField = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.rectHLabel, 2, 0)
        self.settingsView.functionGridLayout.addWidget(self.rectHField, 2, 1)

        self.rectY0Label = QLabel("y0: ", self.settingsView)
        self.rectY0Field = QLineEdit("0", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.rectY0Label, 2, 2)
        self.settingsView.functionGridLayout.addWidget(self.rectY0Field, 2, 3)

        # lin
        self.linY0Label = QLabel("y0: ", self.settingsView)
        self.linY0Field = QLineEdit("0", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.linY0Label, 1, 0)
        self.settingsView.functionGridLayout.addWidget(self.linY0Field, 1, 1)

        self.linY1Label = QLabel("y1: ", self.settingsView)
        self.linY1Field = QLineEdit("1", self.settingsView)
        self.settingsView.functionGridLayout.addWidget(self.linY1Label, 1, 2)
        self.settingsView.functionGridLayout.addWidget(self.linY1Field, 1, 3)

        # preview
        self.previewPlot = SimplePlotter(None, width=5, height=2, dpi=60)
        self.settingsView.functionGridLayout.addWidget(self.previewPlot, 3, 0, 1, 4)

        self.settingsView.previewButton.clicked.connect(self.previewFunction)

    def updateSettings(self):
        self.painter.nsPerFrame = float(self.settingsView.nsPerFrameField.text())
        self.painter.threshold = float(self.settingsView.thresholdField.text())
        self.painter.rendered = False
        self.painter.colorScheme = self.colorScheme
        self.painter.customColor = self.customColor
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def updateFilters(self):
        print("filter update")
        self.painter.labelView.clean()
        # total time filter
        totalTimeActive = self.settingsView.activeTotalTimeCheckbox.isChecked()
        scoreActive = self.settingsView.activeScoreCheckbox.isChecked()
        sortingActive = self.settingsView.activeSortingBox.isChecked()
        filterActive = (totalTimeActive or scoreActive or sortingActive)
        weightActive = self.settingsView.functionActiveCheckbox.isChecked()
        # only filter given range
        rangeFilterActive = self.settingsView.filterRangeCheckbox.isChecked()
        if len(self.contacts) > 0:
            lower = int(self.settingsView.lowerRangeField.text()) - 1
            upper = self.settingsView.upperRangeField.text()
            if upper == "end":
                upper = len(self.contacts[0].scoreArray)
            else:
                upper = int(upper)

            if lower < 0:
                lower = 0
            self.painter.range = [lower, upper]
            self.painter.rangeFilterActive = False
            self.filteredContacts = copy.deepcopy(self.contacts)
            # residue range filter
            resrangeFilter = ResidueRangeFilter("resrange")
            self.filteredContacts = resrangeFilter.filterResiduesByRange(self.filteredContacts, self.settingsView.residARangeField.text(), self.settingsView.residBRangeField.text())
            # aminoacids name filter
            aaFilter = NameFilter("name")
            self.filteredContacts = aaFilter.filterResiduesByName(self.filteredContacts, self.settingsView.residANameField.text(), self.settingsView.residBNameField.text())
            # range filter
            if rangeFilterActive:
                self.painter.rangeFilterActive = True
                frameRangeFilter = FrameFilter("framer")
                self.filteredContacts = frameRangeFilter.extractFrameRange(self.filteredContacts, [lower, upper])
            # weight functions
            if weightActive:
                if self.currentFunctionType == FunctionType.sigmoid:
                    print("sig weight")
                    x0 = float(self.sigX0Field.text())
                    L = float(self.sigLField.text())
                    k = float(self.sigKField.text())
                    y0 = float(self.sigY0Field.text())
                    sig = SigmoidWeightFunction("sig", np.arange(0, len(self.contacts[0].scoreArray), 1), x0, L, k, y0)
                    self.filteredContacts = sig.weightContactFrames(self.filteredContacts)
                elif self.currentFunctionType == FunctionType.rect:
                    x0 = float(self.rectX0Field.text())
                    x1 = float(self.rectX1Field.text())
                    h = float(self.rectHField.text())
                    y0 = float(self.rectY0Field.text())
                    rect = RectangularWeightFunction("rect", np.arange(0, len(self.contacts[0].scoreArray), 1), x0, x1, h, y0)
                    self.filteredContacts = rect.weightContactFrames(self.filteredContacts)
                elif self.currentFunctionType == FunctionType.linear:
                    y0 = float(self.linY0Field.text())
                    y1 = float(self.linY1Field.text())
                    lin = LinearWeightFunction("rect", np.arange(0, len(self.contacts[0].scoreArray), 1), y0, y1)
                    self.filteredContacts = lin.weightContactFrames(self.filteredContacts)
            # other filters
            if filterActive:
                    if totalTimeActive:
                        operator = self.settingsView.compareTotalTimeDropdown.currentText()
                        value = float(self.settingsView.totalTimeField.text())
                        filter = TotalTimeFilter("tottime", operator, value)
                        self.filteredContacts = filter.filterContacts(self.filteredContacts)
                    if scoreActive:
                        operator = self.settingsView.compareScoreDropdown.currentText()
                        value = float(self.settingsView.scoreField.text())
                        filter = ScoreFilter("score", operator, value, self.settingsView.meanDropdown.currentText())
                        self.filteredContacts = filter.filterContacts(self.filteredContacts)
                    if sortingActive:
                        key = self.settingsView.sortingKeyDropdown.currentText()
                        descending = SortingOrder.mapping[self.settingsView.sortingOrderDropdown.currentText()]
                        sorter = Sorting("sorting", key, descending)
                        sorter.setThresholdAndNsPerFrame(float(self.settingsView.thresholdField.text()), float(self.settingsView.nsPerFrameField.text()))
                        self.filteredContacts = sorter.sortContacts(self.filteredContacts)
                    self.painter.contacts = self.filteredContacts
                    self.painter.rendered = False
                    self.painter.update()
                    self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))
                    if len(self.filteredContacts) == 0:
                        self.painter.labelView.clean()
            else:
                #no weight or filters
                self.painter.contacts = self.filteredContacts
                self.painter.rendered = False
                self.painter.update()
                self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    # switch between weight functions
    def showFunctionSettings(self, radiobutton):
        self.currentFunctionType = radiobutton
        if radiobutton == FunctionType.sigmoid:
            self.showHide(False, True, True)
        elif radiobutton == FunctionType.rect:
            self.showHide(True, False, True)
        elif radiobutton == FunctionType.linear:
            self.showHide(True, True, False)
    #hiding and showing of weight function labels and textfields
    def showHide(self, first, second, third):
        self.sigX0Label.setHidden(first)
        self.sigX0Field.setHidden(first)
        self.sigLLabel.setHidden(first)
        self.sigLField.setHidden(first)
        self.sigKLabel.setHidden(first)
        self.sigKField.setHidden(first)
        self.sigY0Label.setHidden(first)
        self.sigY0Field.setHidden(first)
        self.rectX0Label.setHidden(second)
        self.rectX0Field.setHidden(second)
        self.rectX1Label.setHidden(second)
        self.rectX1Field.setHidden(second)
        self.rectHLabel.setHidden(second)
        self.rectHField.setHidden(second)
        self.rectY0Label.setHidden(second)
        self.rectY0Field.setHidden(second)
        self.linY0Label.setHidden(third)
        self.linY0Field.setHidden(third)
        self.linY1Label.setHidden(third)
        self.linY1Field.setHidden(third)

    #display currenct function in preview window
    def previewFunction(self):
        x = []
        y = []
        if self.currentFunctionType == FunctionType.sigmoid:
            x0 = float(self.sigX0Field.text())
            L = float(self.sigLField.text())
            k = float(self.sigKField.text())
            y0 = float(self.sigY0Field.text())
            if len(self.contacts) > 0:
                sig = SigmoidWeightFunction("sig", np.arange(0, len(self.contacts[0].scoreArray), 1), x0, L, k, y0)
                x = np.arange(0, len(self.contacts[0].scoreArray), 1)
                y = sig.previewFunction()
        elif self.currentFunctionType == FunctionType.rect:
            x0 = float(self.rectX0Field.text())
            x1 = float(self.rectX1Field.text())
            h = float(self.rectHField.text())
            y0 = float(self.rectY0Field.text())
            if len(self.contacts) > 0:
                rect = RectangularWeightFunction("rect", np.arange(0, len(self.contacts[0].scoreArray), 1), x0, x1, h, y0)
                x = np.arange(0, len(self.contacts[0].scoreArray), 1)
                y = rect.previewFunction()
        elif self.currentFunctionType == FunctionType.linear:
            y0 = float(self.linY0Field.text())
            y1 = float(self.linY1Field.text())
            if len(self.contacts) > 0:
                lin = LinearWeightFunction("rect", np.arange(0, len(self.contacts[0].scoreArray), 1), y0, y1)
                x = np.arange(0, len(self.contacts[0].scoreArray), 1)
                y = lin.previewFunction()
        sip.delete(self.previewPlot)
        self.previewPlot = SimplePlotter(None, width=5, height=2, dpi=60)
        self.settingsView.functionGridLayout.addWidget(self.previewPlot, 3, 0, 1, 4)
        self.previewPlot.plot(x, y)
        self.previewPlot.update()

    def openPrefs(self):
        self.settingsView.show()

    def showStatistics(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        numberTitleLabel = QLabel("total number of contacts:")
        numberLabel = QLabel(str(len(self.contacts)))

        numberFramesTitleLabel = QLabel("number of frames:")
        numberFramesLabel = QLabel(str(len(self.contacts[0].scoreArray)))

        meanTitleLabel = QLabel("mean contact score:")
        meanLabel = QLabel(str(mean_score_of_contactArray(self.contacts)))

        medianTitleLabel = QLabel("median score:")
        medianLabel = QLabel(str(median_score_of_contactArray(self.contacts)))

        grid.addWidget(numberTitleLabel, 0, 0)
        grid.addWidget(numberLabel, 0, 1)
        grid.addWidget(numberFramesTitleLabel, 1, 0)
        grid.addWidget(numberFramesLabel, 1, 1)
        grid.addWidget(meanTitleLabel, 2, 0)
        grid.addWidget(meanLabel, 2, 1)

        grid.addWidget(medianTitleLabel, 2, 2)
        grid.addWidget(medianLabel, 2, 3)

        allContactPlot = ContactPlotter(None, width=4, height=2, dpi=80)
        allContactPlot.plot_all_contacts_figure(self.contacts)
        grid.addWidget(allContactPlot, 3, 0, 1, 4)
        d.setWindowTitle("Statistics")
        d.resize(600, 450)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def pushOpen(self):
        fnames = QFileDialog.getOpenFileNames(self, "Open file")
        for file in fnames[0]:
            self.file = file
            break
        lines = []
        print(self.file)
        with open(self.file, "r") as f:
            for line in f.readlines():
                lines.append(line)
        self.contacts = makeContactFromLines(lines)
        print("new contacts: " + str(len(self.contacts)))
        self.painter.contacts = self.contacts
        self.painter.range = [0, len(self.contacts[0].scoreArray)]
        # self.painter.update()
        # set max slider value to frame number!
        self.mergeSlider.setMaximum(len(self.contacts[0].scoreArray) / 15)
        self.updateSettings()
        self.updateFilters()

    def pushExport(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        self.exportLabel = QLabel("Export current view: ")

        self.saveButton = QPushButton("Export")
        self.saveButton.setAutoDefault(False)
        self.saveButton.clicked.connect(self.pushSave)

        self.formatBox = QComboBox()
        self.formatBox.addItem("PNG")
        self.formatBox.addItem("SVG")

        grid.addWidget(self.exportLabel, 0, 0)
        grid.addWidget(self.saveButton, 0, 1)

        grid.addWidget(self.formatBox, 2, 0)

        d.setWindowTitle("Export")
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def pushSave(self):
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        print(self.formatBox.currentText())
        if self.formatBox.currentText() == "PNG":
            if len(fileName[0]) > 0:
                print("Saving current view to ", fileName[0])
                currentView = self.painter.grab()
                currentView.save(fileName[0])
        elif self.formatBox.currentText() == "SVG":
            if len(fileName[0]) > 0:
                print("Saving current view to ", fileName[0])
                generator = QSvgGenerator()
                generator.setFileName(fileName[0])
                generator.setSize(self.painter.size())
                generator.setViewBox(self.painter.rect())
                self.painter.renderContact(generator)
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def mergeValueChanged(self):
        self.painter.merge = self.mergeSlider.value()
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def alphaValueChanged(self):
        self.painter.alphaFactor = self.alphaSlider.value()
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def showColorPicker(self):
        col = QColorDialog.getColor()
        self.customColor = col
        if col.isValid():
            self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }"
                                   % self.customColor.name())

    def updateColors(self):
        if self.settingsView.bbscScoreRadioButton.isChecked():
            self.colorScheme = ColorScheme.bbsc
        elif self.settingsView.customColorRadioButton.isChecked():
            self.colorScheme = ColorScheme.custom
        self.updateSettings()
        self.updateFilters()

    def visualize(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        label = QLabel("Split selections for each contact")
        self.splitVisCheckbox = QCheckBox()
        grid.addWidget(label, 0, 0)
        grid.addWidget(self.splitVisCheckbox, 0, 1)

        self.additionalText1 = QLineEdit()
        self.additionalText2 = QLineEdit()
        additionalTextLabel1 = QLabel("Additional seltext for selection 1")
        additionalTextLabel2 = QLabel("Additional seltext for selection 2")

        button = QPushButton("Create tcl script")
        button.clicked.connect(self.createTclScriptVis)
        grid.addWidget(button,3,0,1,2)

        grid.addWidget(additionalTextLabel1, 1, 0)
        grid.addWidget(additionalTextLabel2, 2, 0)
        grid.addWidget(self.additionalText1, 1, 1)
        grid.addWidget(self.additionalText2, 2, 1)

        d.setWindowTitle("Visualize in VMD")
        d.resize(200, 300)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def createTclScriptVis(self):
        f = open('vis.tcl', 'w')
        f.write('mol new %s \n' % self.psf)
        f.write('mol addfile %s \n' % self.dcd)
        f.write('mol delrep 0 top \n')
        f.write('mol representation NewCartoon \n')
        f.write('mol Color ColorID 3 \n')
        f.write('mol selection {all} \n')
        f.write('mol addrep top \n')

        if self.splitVisCheckbox.isChecked():
            for cont in self.filteredContacts:
                currentSel1 = []
                index = 0
                for item in cont.key1:
                    if item != "none":
                        currentSel1.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                    index += 1
                currentSel1String = "and".join(currentSel1)
                currentSel2 = []
                index = 0
                for item in cont.key2:
                    if item != "none":
                        currentSel2.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                    index += 1
                currentSel2String = " and ".join(currentSel2)
                add1 = ("" if self.additionalText1.text()=="" else (" and " + self.additionalText1.text()))
                add2 = ("" if self.additionalText2.text() == "" else (" and " + self.additionalText2.text()))
                sel = "("+currentSel1String + add1 + ") or (" + currentSel2String + add2 + ")"
                f.write('mol representation Licorice \n')
                f.write('mol Color Name \n')
                f.write('mol selection {%s} \n'%sel)
                f.write('mol addrep top \n')
        else:
            total = []
            for cont in self.filteredContacts:
                currentSel1 = []
                index = 0
                for item in cont.key1:
                    if item != "none":
                        currentSel1.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                    index += 1
                currentSel1String = "and".join(currentSel1)
                currentSel2 = []
                index = 0
                for item in cont.key2:
                    if item != "none":
                        currentSel2.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                    index += 1
                currentSel2String = " and ".join(currentSel2)
                add1 = ("" if self.additionalText1.text() == "" else (" and " + self.additionalText1.text()))
                add2 = ("" if self.additionalText2.text() == "" else (" and " + self.additionalText2.text()))
                sel = "(" + currentSel1String + add1 + ") or (" + currentSel2String + add2 + ")"
                total.append(sel)
            seltext = " or ".join(total)
            f.write('mol representation Licorice \n')
            f.write('mol Color Name \n')
            f.write('mol selection {%s} \n' % seltext)
            f.write('mol addrep top \n')
        f.close()



class SettingsTabWidget(QTabWidget, Ui_settingsWindowWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QTabWidget, self).__init__(parent)
        self.setupUi(self)


class ColorScheme:
    custom, bbsc = range(2)