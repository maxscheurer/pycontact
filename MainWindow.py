'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
from __future__ import print_function
import multiprocessing
import warnings
import time
import itertools
import pickle
import copy

from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QProgressBar
from PyQt5.QtGui import QPaintEvent
from PyQt5.Qt import Qt
import numpy as np
from numpy import linalg as la
from matplotlib import cm

import gui
from multi_accumulation import *
from biochemistry import vdwRadius
from SasaWidgets import *
from Canvas import *
from Dialogues import FileLoaderDialog, AnalysisDialog
from ExportTabWidget import ExportTabWidget
from Plotters import *
from ContactAnalyzer import *
from ErrorBox import ErrorBox
from LogPool import *
from aroundPatch import AroundSelection

multiprocessing.log_to_stderr()
np.set_printoptions(threshold=np.inf)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")


class MainWindow(QMainWindow, gui.Ui_MainWindow):
    def __init__(self, parent=None):
        self.config = None
        self.analysis = None
        self.maps = None
        super(MainWindow, self).__init__(parent)
        self.contacts = []
        self.filteredContacts = []
        self.setupUi(self)

        self.setWindowTitle("pyContact")
        self.mergeSlider.valueChanged.connect(self.mergeValueChanged)

        # painter contains both labels and frame boxes for drawing
        self.painter = Canvas()
        self.scrollArea.setWidget(self.painter)
        # deprecated
        # self.actionOpen.triggered.connect(self.pushOpen)
        self.actionExport.triggered.connect(self.pushExport)
        self.actionLoad_Data.triggered.connect(self.loadDataPushed)
        self.actionExport_Session.triggered.connect(self.exportSession)
        self.actionImport_Session.triggered.connect(self.importSession)
        self.actionShow_Info.triggered.connect(self.showDeveloperInfo)
        # settings and filters
        self.settingsView = SettingsTabWidget()
        self.settingsView.applySettingsButton.clicked.connect(self.updateSettings)
        self.settingsView.applyFilterButton.clicked.connect(self.updateFilters)

# deprecated
        #alpha slider for color
        # self.alphaSlider.setValue(50)
        # self.alphaSlider.valueChanged.connect(self.alphaValueChanged)

        self.statisticsButton.clicked.connect(self.showStatistics)

        self.openPreferencesButton.clicked.connect(self.openPrefs)

        #color picker
        self.settingsView.pickColorButton.clicked.connect(self.showColorPicker)
        self.customColor = QColor(230, 50, 0)
        self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" % self.customColor.name())

        #analysis button
        self.analysisButton.clicked.connect(self.analyzeDataPushed)


        #contact area button
        self.contactAreaButton.clicked.connect(self.showContactAreaView)

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

        self.actionDefault.triggered.connect(self.loadDefault)

        self.progressWidget = ProgessWidget("Progress")

        self.exportWidget = ExportTabWidget()

        self.sasaView = SasaWidget()

        self.updateSettings()
        self.updateFilters()

    def importSession(self):
        fnames = QFileDialog.getOpenFileNames(self, "Open file")
        importfile = ""
        for file in fnames[0]:
            importfile = file
            break
        importDict = pickle.load(open(importfile, "rb"))
        self.contacts = importDict["contacts"]
        arguments = importDict["analyzer"][0:-1]
        trajArgs = importDict["trajectory"]
        self.analysis = Analyzer(*arguments)
        self.analysis.contactResults = importDict["analyzer"][-1]
        self.analysis.setTrajectoryData(*trajArgs)
        self.updateSettings()
        self.updateFilters()

    def exportSession(self):
        fileName = QFileDialog.getSaveFileName(self, 'Export file')
        filestring = fileName[0]
        if filestring == "":
            return
        if self.contacts is not None and self.analysis is not None:
            analyzerArgs = [self.analysis.psf, self.analysis.dcd, self.analysis.cutoff, self.analysis.hbondcutoff, self.analysis.hbondcutangle, self.analysis.sel1text, self.analysis.sel2text,self.analysis.contactResults]
            trajArgs = self.analysis.getTrajectoryData()
            exportDict = {"contacts":self.contacts,"analyzer":analyzerArgs,"trajectory":trajArgs}
            pickle.dump(exportDict, open(filestring, "wb"))
        else:
            box = ErrorBox("No data to export.")
            box.exec_()
            return

    def loadDefault(self):
        # importDict = pickle.load(open("defaultsession", "rb"))
        importDict = pickle.load(open("defaultsession", "rb"))
        self.contacts = importDict["contacts"]
        arguments = importDict["analyzer"][0:-1]
        trajArgs = importDict["trajectory"]
        self.analysis = Analyzer(*arguments)
        self.analysis.contactResults = importDict["analyzer"][-1]
        self.analysis.setTrajectoryData(*trajArgs)
        self.updateSettings()
        self.updateFilters()

    def loadData_parallel(self, nprocs):
        from multi_trajectory import run_load_parallel
        # run_load_parallel(nproc, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text)
        return run_load_parallel(nprocs,self.config.psf,self.config.dcd, self.config.cutoff,self.config.hbondcutoff,self.config.hbondcutangle,self.config.sel1text,self.config.sel2text)

    def loadDataPushed(self):
        self.config,result = FileLoaderDialog.getConfig()
        if result == 1:
            attrs = vars(self.config)
            nproc = int(self.settingsView.coreBox.value())
            # do not allow multiprocessing unless the trajectory has enough frames
            if nproc == 1:
                parallel = 0
            else:
                parallel = 1
            print(', '.join("%s: %s" % item for item in attrs.items()))
            self.analysis = Analyzer(self.config.psf,self.config.dcd, self.config.cutoff,self.config.hbondcutoff,self.config.hbondcutangle,self.config.sel1text,self.config.sel2text)
            # self.connect(self.analysis, QtCore.SIGNAL('taskUpdated'),self.handleTaskUpdated)
            # self.connect(self.analysis, QtCore.SIGNAL('frameNumberSet'),self.setFrameNumber)
            # self.progressWidget.show()
            if parallel:
                self.analysis.contactResults,self.analysis.resname_array,self.analysis.resid_array,self.analysis.name_array,self.analysis.type_array,self.analysis.segids,self.analysis.backbone,self.analysis.sel1text,self.analysis.sel2text = self.loadData_parallel(nproc)
            else:
                self.analysis.runFrameScan()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Data loaded: %f frames scanned." % len(self.analysis.contactResults))
            msg.setInformativeText("")
            msg.setWindowTitle("Data loaded")
            msg.setDetailedText("Now click on Analysis to proceed")
            msg.exec_()

    # progress of loading trajectory
    def handleTaskUpdated(self):
    	print(self.analysis.currentFrame)
    	self.progressWidget.setValue(self.analysis.currentFrame)

    # progress of loading trajectory
    def setFrameNumber(self):
        self.progressWidget.setMax(self.analysis.totalFrameNumber)

    def analyzeParallel(self, map1, map2):
        nproc = int(self.settingsView.coreBox.value())
        start = time.time()
        trajData = self.analysis.getTrajectoryData()
        contResults = self.analysis.contactResults
        results = []
        rank = 0
        manager = multiprocessing.Manager()
        d = manager.list(trajData)
        all_chunk = chunks(contResults, nproc)
        pool = LoggingPool(nproc)
        print("Running on %d cores" % nproc)
        for c in all_chunk:
            results.append(pool.apply_async(loop_frame, args=(c, map1, map2, d)))
            rank += 1
        # TODO: might be important, but without, it's faster and until now, provides the same results
        pool.close()
        pool.join()
        stop = time.time()
        print("time: ", str(stop-start), rank)
        print(str(len(c)), rank)
        allkeys = []
        frame_contacts_accumulated = []
        print(len(results))
        for res in results:
            rn = res.get()
            allkeys.extend(rn[0])
            frame_contacts_accumulated.extend(rn[1])
        accumulatedContactsDict = {}
        #   start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if key not in frame_dict:  # puts empty score TempContactAccumulate in dict
                    key1, key2 = makeKeyArraysFromKey(key)
                    emptyCont = TempContactAccumulate(key1, key2)
                    emptyCont.fscore = 0
                    frame_dict[key] = emptyCont
                accumulatedContactsDict[key].append(frame_dict[key])
        finalAccumulatedContacts = []  # list of AccumulatedContacts
        for key in accumulatedContactsDict:
            key1, key2 = makeKeyArraysFromKey(key)
            acc = AccumulatedContact(key1, key2)
            for tempContact in accumulatedContactsDict[key]:
                acc.addScore(tempContact.fscore)
                acc.addContributingAtoms(tempContact.contributingAtomContacts)
                acc.bb1 += tempContact.bb1score
                acc.bb2 += tempContact.bb2score
                acc.sc1 += tempContact.sc1score
                acc.sc2 += tempContact.sc2score
            finalAccumulatedContacts.append(acc)
        # stop = time.time()
        # print(stop - start)
        glob_stop = time.time()
        print(glob_stop - start)
        return finalAccumulatedContacts

    def analyzeDataPushed(self):
        self.maps, result = AnalysisDialog.getMapping()
        if result == 1:
            map1 = self.maps[0]
            map2 = self.maps[1]
            nproc = int(self.settingsView.coreBox.value())
            if nproc == 1:
                parallel = 0
            else:
                parallel = 1

            if parallel:
                self.contacts = self.analyzeParallel(map1, map2)
            else:
                self.contacts = self.analysis.runContactAnalysis(map1, map2)
            # for cont in self.contacts:
                # cont.setScores()
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
        self.painter.showHbondScores = False
        # total time filter
        totalTimeActive = self.settingsView.activeTotalTimeCheckbox.isChecked()
        scoreActive = self.settingsView.activeScoreCheckbox.isChecked()
        sortingActive = self.settingsView.activeSortingBox.isChecked()
        onlyActive = self.settingsView.onlyBoxActiveCheckbox.isChecked()
        filterActive = (totalTimeActive or scoreActive or sortingActive or onlyActive)
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
            range_filter = RangeFilter("resrange")
            self.filteredContacts = range_filter.filterByRange(self.filteredContacts, self.settingsView.residARangeField.text(), self.settingsView.residBRangeField.text(),AccumulationMapIndex.resid)

            self.filteredContacts = range_filter.filterByRange(self.filteredContacts, self.settingsView.atomAIndexField.text(), self.settingsView.atomBIndexField.text(),AccumulationMapIndex.index)

            # aminoacids name filter
            name_filter = NameFilter("name")
            self.filteredContacts = name_filter.filterContactsByName(self.filteredContacts, self.settingsView.residANameField.text(), self.settingsView.residBNameField.text(),AccumulationMapIndex.resname)

            self.filteredContacts = name_filter.filterContactsByName(self.filteredContacts, self.settingsView.atomANameField.text(), self.settingsView.atomBNameField.text(),AccumulationMapIndex.name)
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
                    if onlyActive:
                        key = self.settingsView.selectOnlyToolbox.currentText()
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

    # switch between weight functions
    def showFunctionSettings(self, radiobutton):
        self.currentFunctionType = radiobutton
        if radiobutton == FunctionType.sigmoid:
            self.showHide(False, True, True)
        elif radiobutton == FunctionType.rect:
            self.showHide(True, False, True)
        elif radiobutton == FunctionType.linear:
            self.showHide(True, True, False)

    # hiding and showing of weight function labels and textfields
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

    # display currenct function in preview window
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
        if len(self.contacts) == 0 or self.contacts == None:
            box = ErrorBox("No data loaded!")
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
        copyright = QLabel("Version 0.1.1a, January 2017")

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
            self.exportWidget.setMaps(self.maps[0],self.maps[1])
            self.exportWidget.setMapLabels(self.analysis.sel1text,self.analysis.sel2text)
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

    def mergeValueChanged(self):
        self.painter.merge = self.mergeSlider.value()
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

# deprecated
    def alphaValueChanged(self):
        self.painter.alphaFactor = 50
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def showColorPicker(self):
        col = QColorDialog.getColor()
        self.customColor = col
        if col.isValid():
            self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" % self.customColor.name())

    def updateColors(self):
        if self.settingsView.bbscScoreRadioButton.isChecked():
            self.colorScheme = ColorScheme.bbsc
        elif self.settingsView.customColorRadioButton.isChecked():
            self.colorScheme = ColorScheme.custom
        self.updateSettings()
        self.updateFilters()

    def showContactAreaView(self):
        self.sasaView.show()


class SettingsTabWidget(QTabWidget, Ui_settingsWindowWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QTabWidget, self).__init__(parent)
        self.setupUi(self)


class ColorScheme:
    custom, bbsc = range(2)
