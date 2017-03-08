from __future__ import print_function
import multiprocessing
import warnings
import copy

import PyQt5.QtCore as QtCore
from PyQt5.QtCore import QRect, pyqtSlot, QObject
from PyQt5.QtWidgets import (QMainWindow, QTabWidget, QLabel, QDialog,
                             QApplication, QGridLayout, QFileDialog, QColorDialog)
from PyQt5.QtGui import QPaintEvent
from PyQt5.Qt import Qt, QColor
from PyQt5.QtGui import QIntValidator
import numpy as np

from . import MainQtGui
from ..core.multi_accumulation import *
from ..core.Biochemistry import vdwRadius
from SasaWidgets import SasaWidget
from Canvas import Canvas
from Dialogues import FileLoaderDialog, AnalysisDialog
from ExportTabWidget import ExportTabWidget
from Plotters import *
from ..core.ContactAnalyzer import *
from ErrorBox import ErrorBox
from ErrorMessages import ErrorMessages
from ..core.LogPool import *
from ..core.aroundPatch import AroundSelection
import settings
from ..exampleData.datafiles import DEFAULTSESSION
from VMDControlPanel import VMDControlPanel
from ..core.DataHandler import DataHandler

multiprocessing.log_to_stderr()
np.set_printoptions(threshold=np.inf)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")


class MainWindow(QMainWindow, MainQtGui.Ui_MainWindow, QObject):
    """PyContact Application Main Window with timeline"""

    def closeEvent(self, event):
        """ Closing application when Exit on MainWindow is clicked"""
        print("Closing Application")
        event.accept()
        QApplication.quit()

    def __init__(self, parent=None):
        self.config = None
        self.analysis = None
        self.maps = None
        super(MainWindow, self).__init__(parent)
        self.contacts = []
        self.filteredContacts = []
        self.setupUi(self)

        self.setWindowTitle("pyContact")

        # painter contains both labels and frame boxes for drawing
        self.painter = Canvas()
        self.scrollArea.setWidget(self.painter)
        self.actionExportData.triggered.connect(self.pushExport)
        self.actionLoad_Data.triggered.connect(self.loadDataPushed)
        self.actionExport_Session.triggered.connect(self.exportSession)
        self.actionImport_Session.triggered.connect(self.importSession)
        self.actionShow_Info.triggered.connect(self.showDeveloperInfo)
        # settings and filters
        self.settingsView = SettingsTabWidget()
        self.settingsView.applySettingsButton.clicked.connect(self.updateSettings)
        self.applyFilterButton.clicked.connect(self.updateFilters)
        # statistics
        self.statisticsButton.clicked.connect(self.showStatistics)
        # color picker
        self.settingsView.pickColorButton.clicked.connect(self.showColorPicker)
        self.customColor = QColor(230, 50, 0)
        self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" % self.customColor.name())

        # frames stride
        posIntValidator = QIntValidator()
        posIntValidator.setBottom(1)
        self.frameStrideField.setValidator(posIntValidator)

        # analysis button
        self.analysisButton.clicked.connect(self.analyzeDataPushed)
        # contact area button
        self.actionContact_Area_Calculations.triggered.connect(self.showContactAreaView)
        # preferences
        self.actionPreferences.triggered.connect(self.openPrefs)
        # apply color button
        self.settingsView.applyColorButton.clicked.connect(self.updateColors)
        self.colorScheme = ColorScheme.bbsc

        self.actionDefault.triggered.connect(self.loadDefault)

        self.currentSelection1 = "-"
        self.currentSelection2 = "-"

        # setup of extra widgets
        self.exportWidget = ExportTabWidget()
        self.sasaView = SasaWidget()

        self.updateSettings()
        self.updateFilters()

        self.analysis_state = False

        self.vismode = False
        self.visModeButton.setCheckable(True)
        self.visModeButton.setChecked(False)
        self.visModeButton.clicked.connect(self.switchedToVisMode)

        self.vmdpanel = VMDControlPanel()
        self.actionVMD_Remote_Control.triggered.connect(self.showVMDControlPanel)

        self.painter.clickedRowSignal.connect(self.updateVMDSelections)
        self.painter.clickedColumnSignal.connect(self.updateVMDFrame)

    def showVMDControlPanel(self):
        self.vmdpanel.show()

    def showContactAreaView(self):
        self.sasaView.show()

    def switchedToVisMode(self):
        if self.visModeButton.isChecked():
            self.vismode = True
        else:
            self.vismode = False
        self.painter.switchToVisMode(self.vismode)
        self.updateSettings()
        self.updateFilters()

    @pyqtSlot()
    def updateVMDSelections(self):
        if self.vmdpanel.connected:
            self.vmdpanel.updateSelections(self.analysis.sel1text, self.analysis.sel2text,
                                           [self.contacts[self.painter.globalClickedRow]])

    @pyqtSlot()
    def updateVMDFrame(self):
        if self.vmdpanel.connected:
            self.vmdpanel.gotoVMDFrame(self.painter.clickedColumn)

    def updateSelectionLabels(self, sel1, sel2):
        self.currentSelection1 = sel1
        self.currentSelection2 = sel2
        self.selection1label.setText(sel1)
        self.selection2label.setText(sel2)

    def importSession(self):
        fnames = QFileDialog.getOpenFileNames(self, "Open file")
        importfile = ""
        for f in fnames[0]:
            importfile = f
            break
        if importfile == "" or len(fnames) == 0:
            return

        self.contacts, arguments, trajArgs, self.maps, contactResults = DataHandler.importSessionFromFile(importfile)
        self.analysis = Analyzer(*arguments)
        self.analysis.contactResults = contactResults
        self.analysis.setTrajectoryData(*trajArgs)
        self.updateSelectionLabels(arguments[5], arguments[6])
        self.updateSettings()
        self.updateFilters()

    def exportSession(self):
        fileName = QFileDialog.getSaveFileName(self, 'Export file')
        filestring = fileName[0]
        if filestring == "":
            return
        if self.contacts is not None and self.analysis is not None:
            self.setInfoLabel("Exporting current session...")
            analyzerArgs = [self.analysis.psf, self.analysis.dcd, self.analysis.cutoff, self.analysis.hbondcutoff,
                            self.analysis.hbondcutangle, self.analysis.sel1text, self.analysis.sel2text,
                            self.analysis.contactResults]
            trajArgs = self.analysis.getTrajectoryData()
            exportDict = {"contacts": self.contacts, "analyzer": analyzerArgs, "trajectory": trajArgs,
                          "maps": [self.analysis.lastMap1, self.analysis.lastMap2]}
            pickle.dump(exportDict, open(filestring, "wb"))
            self.cleanInfoLabel()
        else:
            box = ErrorBox(ErrorMessages.NOEXPDATA)
            box.exec_()
            return

    def loadDefault(self):
        self.contacts, arguments, trajArgs, self.maps, contactResults = \
            DataHandler.importSessionFromFile(DEFAULTSESSION)
        self.analysis = Analyzer(*arguments)
        self.analysis.contactResults = contactResults
        self.analysis.setTrajectoryData(*trajArgs)
        self.updateSelectionLabels(arguments[5], arguments[6])
        self.updateSettings()
        self.updateFilters()

    def loadDataPushed(self):
        self.config, result = FileLoaderDialog.getConfig()
        if result == 1:
            self.setInfoLabel("Loading trajectory and running atomic contact analysis...")
            nproc = int(self.settingsView.coreBox.value())
            self.analysis = Analyzer(self.config.psf, self.config.dcd, self.config.cutoff, self.config.hbondcutoff,
                                     self.config.hbondcutangle, self.config.sel1text, self.config.sel2text)
            QApplication.processEvents()
            self.analysis.runFrameScan(nproc)
            self.setInfoLabel("%d frames loaded." % len(self.analysis.contactResults))
            self.updateSelectionLabels(self.config.sel1text, self.config.sel2text)

    # progress of loading trajectory
    # def handleTaskUpdated(self):
    #     self.progressBar.setValue(self.analysis.currentFrame)

    # progress of loading trajectory
    # def setFrameNumber(self):
    #     self.progressBar.setMax(self.analysis.totalFrameNumber)

# deprecated?
    # @pyqtSlot()
    # def updateAnalyzedFrames(self):
    #     QApplication.processEvents()
    #     self.progressBar.setValue(100* float(self.value) / float(self.totalFramesToProcess))
    #     self.value += 1

    # def analysisEventListener(self):
    #     while self.analysis_state:
    #         QApplication.processEvents()
    #         progress = 0
    #         for each in analysisProgressDict.keys():
    #             progress += analysisProgressDict[each]
    #             # sasaProgressDict[each] = 0
    #         progress = float(progress) / float(self.totalFramesToProcess) * 100
    #         # if (101 - self.sasaProgressBar.value()) < progress:
    #         #     self.sasaProgressBar.update_bar(101 - self.sasaProgressBar.value())
    #         if progress > 0:
    #             # print(progress)
    #             self.progressBar.setValue(progress)
    #
    #         if int(progress) == 100:
    #             # print("finished")
    #             for each in analysisProgressDict.keys():
    #                 analysisProgressDict[each]=0
    #             progress = 0
    #             self.progressBar.setValue(0)
    #             self.analysis_state = False

    def setInfoLabel(self, txt):
        self.statusLabel.setText(txt)

    def cleanInfoLabel(self):
        self.setInfoLabel("-")

    def analyzeDataPushed(self):
        if self.analysis is None:
            box = ErrorBox(ErrorMessages.NODATA_PROMPTLOAD)
            box.exec_()
            return

        self.maps, result = AnalysisDialog.getMapping()
        if result == 1:
            self.setInfoLabel("Analyzing contacts...")
            map1 = self.maps[0]
            map2 = self.maps[1]
            nproc = int(self.settingsView.coreBox.value())
            self.contacts = self.analysis.runContactAnalysis(map1, map2, nproc)
            self.progressBar.setValue(0)
            self.setInfoLabel("Updating timeline...")
            QApplication.processEvents()
            self.updateSettings()
            self.updateFilters()
            self.cleanInfoLabel()

    def updateSettings(self):
        self.painter.nsPerFrame = float(self.settingsView.nsPerFrameField.text())
        self.painter.threshold = float(self.settingsView.thresholdField.text())
        self.painter.rendered = False
        self.painter.colorScheme = self.colorScheme
        self.painter.customColor = self.customColor
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def updateFilters(self):
        stride = int(self.frameStrideField.text())
        if stride < 1:
            stride = 1
            QApplication.processEvents()
            self.frameStrideField.setText(str(stride))
        print("stride: ", stride)
        self.painter.merge = stride
        self.painter.labelView.clean()
        self.painter.showHbondScores = False
        # total time filter
        totalTimeActive = self.activeTotalTimeCheckbox.isChecked()
        scoreActive = self.activeScoreCheckbox.isChecked()
        sortingActive = self.activeSortingBox.isChecked()
        onlyActive = self.onlyBoxActiveCheckbox.isChecked()
        filterActive = (totalTimeActive or scoreActive or sortingActive or onlyActive)
        weightActive = False
        # only filter given range
        rangeFilterActive = self.filterRangeCheckbox.isChecked()
        if len(self.contacts) > 0:
            lower = int(self.lowerRangeField.text()) - 1
            upper = self.upperRangeField.text()
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
            range_filter = RangeFilter("resrange")
            self.filteredContacts = range_filter.filterByRange(self.filteredContacts, self.residARangeField.text(),
                                                               self.residBRangeField.text(), AccumulationMapIndex.resid)

            self.filteredContacts = range_filter.filterByRange(self.filteredContacts, self.atomAIndexField.text(),
                                                               self.atomBIndexField.text(), AccumulationMapIndex.index)

            # aminoacids name filter
            name_filter = NameFilter("name")
            self.filteredContacts = name_filter.filterContactsByName(self.filteredContacts, self.residANameField.text(),
                                                                     self.residBNameField.text(),
                                                                     AccumulationMapIndex.resname)

            self.filteredContacts = name_filter.filterContactsByName(self.filteredContacts, self.atomANameField.text(),
                                                                     self.atomBNameField.text(),
                                                                     AccumulationMapIndex.name)
            # range filter
            if rangeFilterActive:
                self.painter.rangeFilterActive = True
                frameRangeFilter = FrameFilter("framer")
                self.filteredContacts = frameRangeFilter.extractFrameRange(self.filteredContacts, [lower, upper])
            for c in self.filteredContacts:
                c.setScores()
                c.setContactType()
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
                        operator = self.compareTotalTimeDropdown.currentText()
                        value = float(self.totalTimeField.text())
                        filter = TotalTimeFilter("tottime", operator, value)
                        self.filteredContacts = filter.filterContacts(self.filteredContacts)
                    if scoreActive:
                        operator = self.compareScoreDropdown.currentText()
                        value = float(self.scoreField.text())
                        filter = ScoreFilter("score", operator, value, self.meanDropdown.currentText())
                        self.filteredContacts = filter.filterContacts(self.filteredContacts)
                    if sortingActive:
                        key = self.sortingKeyDropdown.currentText()
                        descending = SortingOrder.mapping[self.sortingOrderDropdown.currentText()]
                        sorter = Sorting("sorting", key, descending)
                        sorter.setThresholdAndNsPerFrame(float(self.thresholdField.text()),
                                                         float(self.nsPerFrameField.text()))
                        self.filteredContacts = sorter.sortContacts(self.filteredContacts)
                    if onlyActive:
                        key = self.selectOnlyToolbox.currentText()
                        only = OnlyFilter("only", key, 0)
                        self.filteredContacts = only.filterContacts(self.filteredContacts)
                        if key == "hbonds":
                            self.painter.showHbondScores = True
                    self.painter.contacts = self.filteredContacts
                    self.painter.rendered = False
                    self.painter.update()
                    self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))
                    if len(self.filteredContacts) == 0:
                        self.painter.labelView.clean()
            else:
                # no weight or filters
                self.painter.showHbondScores = False
                self.painter.contacts = self.filteredContacts
                self.painter.rendered = False
                self.painter.update()
                self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

        # Update data for export
        self.exportWidget.setContacts(self.filteredContacts)
        if self.maps is not None:
            self.exportWidget.setMaps(self.maps[0], self.maps[1])
            self.exportWidget.setMapLabels(self.analysis.sel1text, self.analysis.sel2text)
        self.exportWidget.setThresholdAndNsPerFrame(self.painter.threshold, self.painter.nsPerFrame)

    def openPrefs(self):
        self.settingsView.show()

    def showStatistics(self):
        if len(self.contacts) == 0 or self.contacts is None:
            box = ErrorBox(ErrorMessages.NOSCORES_PROMPTANALYSIS)
            box.exec_()
            return
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

    def showDeveloperInfo(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        info = QLabel("Developers: Maximilian Scheurer and Peter Rodenkirch")
        info2 = QLabel("Departments: TCBG, University of Illinois at Urbana-Champaign; BZH Heidelberg University")
        mail = QLabel("Contact: mscheurer@ks.uiuc.edu, rodenkirch@stud.uni-heidelberg.de")
        copyright = QLabel("Version 0.1.1a, March 2017")

        grid.addWidget(info, 0, 0)
        grid.addWidget(info2, 1, 0)
        grid.addWidget(mail, 2, 0)
        grid.addWidget(copyright, 3, 0)

        d.setWindowTitle("Developer Info")
        d.resize(150, 80)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def pushExport(self):
        self.exportWidget.valueUpdated.connect(self.handleExportUpdate)
        self.exportWidget.setContacts(self.filteredContacts)
        if self.maps is not None:
            self.exportWidget.setMaps(self.maps[0], self.maps[1])
            self.exportWidget.setMapLabels(self.analysis.sel1text, self.analysis.sel2text)
        self.exportWidget.setThresholdAndNsPerFrame(self.painter.threshold, self.painter.nsPerFrame)
        self.exportWidget.show()

    @QtCore.Slot(str, str)
    def handleExportUpdate(self, fileName, fileType):
        if fileType == "PNG":
            if len(fileName) > 0:
                print("Saving current view to ", fileName)
                currentView = self.painter.grab()
                currentView.save(fileName)
        elif fileType == "SVG":
            if len(fileName) > 0:
                print("Saving current view to ", fileName)
                generator = QSvgGenerator()
                generator.setFileName(fileName)
                generator.setSize(self.painter.size())
                generator.setViewBox(self.painter.rect())
                self.painter.renderContact(generator)
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

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

    def showColorPicker(self):
        col = QColorDialog.getColor()
        self.customColor = col
        if col.isValid():
            self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" %
                                                            self.customColor.name())

    def updateColors(self):
        if self.settingsView.bbscScoreRadioButton.isChecked():
            self.colorScheme = ColorScheme.bbsc
        elif self.settingsView.customColorRadioButton.isChecked():
            self.colorScheme = ColorScheme.custom
        self.updateSettings()
        self.updateFilters()


class SettingsTabWidget(QTabWidget, settings.Ui_settingsWindowWidget):
    def __init__(self, parent=None):
        super(QTabWidget, self).__init__(parent)
        self.setupUi(self)


class ColorScheme:
    custom, bbsc = range(2)
