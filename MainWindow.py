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
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QProgressBar
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore");
    import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import pickle
from matplotlib.mlab import bivariate_normal
from mpl_toolkits.mplot3d import Axes3D

import multiprocessing
from multi_accumulation import *
from sasa_gui import *


sasaProgressManager = multiprocessing.Manager()
sasaProgressDict = sasaProgressManager.dict()

class MainWindow(QMainWindow, gui.Ui_MainWindow):
    def __init__(self, parent=None):
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

        #alpha slider for color
        self.alphaSlider.setValue(50)
        self.alphaSlider.valueChanged.connect(self.alphaValueChanged)

        self.statisticsButton.clicked.connect(self.showStatistics)

        self.openPreferencesButton.clicked.connect(self.openPrefs)

        #color picker
        self.settingsView.pickColorButton.clicked.connect(self.showColorPicker)
        self.customColor = QColor(230, 50, 0)
        self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" % self.customColor.name())

        #analysis button
        self.analysisButton.clicked.connect(self.analyzeDataPushed)


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

        self.actionDefault.triggered.connect(self.loadDefault)

        self.progressWidget = ProgessWidget("Progress")

        self.sasaView = SasaWidget()
        self.sasaView.show()
        # map1 = [0, 0, 0, 1, 1, 0]
        # map2 = [0, 0, 0, 1, 1, 0]
        # self.map1 = map1
        # self.map2 = map2
        # self.psf = "rpn11_ubq_interface-ionized.psf"
        # self.dcd = "short.dcd"
        # analysis.runFrameScan()
        # self.contacts= analysis.runContactAnalysis(map1, map2)

        # map1 = [0, 0, 0, 1, 1, 0]
        # map2 = [0, 0, 0, 1, 1, 0]
        # contactResults = analyze_psf_dcd("membraneAligned-helixZ90X0-ionized.psf", "arf_hmmm_100ns_100f_comet.dcd", 5.0, 2.5, 120, "segid PROT","segid MEMB")
        # self.contacts= analyze_contactResultsWithMaps(contactResults, map1, map2)

        # testing shit
        # maxresids1 = []
        # maxresids2 = []
        # for cont in self.contacts:
        #     cont.determineBackboneSidechainType()
        #     maxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
        #     maxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))

        # # Generate some test data
        # x = np.arange(1,np.max(maxresids1)+2)
        # y = np.arange(1,np.max(maxresids2)+2)
        # data = np.zeros((len(x), len(y)))
        # for cont in self.contacts:
        #     r1 = int(cont.key1[AccumulationMapIndex.resid])
        #     r2 = int(cont.key2[AccumulationMapIndex.resid])
        #     hbonds = cont.hbondFramesScan()
        #     count = np.count_nonzero(hbonds)
        #     if count > 0:
        #         print cont.title + " contains " + str(count) + " hbonds in total"
        #     # data[r1,r2] = cont.total_time(1, 0)
        #     data[r1, r2] = cont.mean_score()
        # fig = plt.figure()
        # ax = fig.add_subplot(111)

        # cax = ax.matshow(data, cmap=cm.coolwarm, aspect='equal')
        # fig.tight_layout()
        # fig.colorbar(cax)
        # plt.show()

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
        exportfile = ""
        filestring = fileName[0]
        if self.contacts is not None and self.analysis is not None:
            analyzerArgs = [self.analysis.psf, self.analysis.dcd, self.analysis.cutoff, self.analysis.hbondcutoff, self.analysis.hbondcutangle, self.analysis.sel1text, self.analysis.sel2text,self.analysis.contactResults]
            trajArgs = self.analysis.getTrajectoryData()
            exportDict = {"contacts":self.contacts,"analyzer":analyzerArgs,"trajectory":trajArgs}
            pickle.dump(exportDict, open(filestring, "wb"))

    def loadDefault(self):
        # importDict = pickle.load(open("defaultsession", "rb"))
        importDict = pickle.load(open("arfsession", "rb"))
        self.contacts = importDict["contacts"]
        arguments = importDict["analyzer"][0:-1]
        trajArgs = importDict["trajectory"]
        self.analysis = Analyzer(*arguments)
        self.analysis.contactResults = importDict["analyzer"][-1]
        self.analysis.setTrajectoryData(*trajArgs)
        self.updateSettings()
        self.updateFilters()

    def loadData_parallel(self):
        from multi_trajectory import run_load_parallel



    def loadDataPushed(self):
        parallel = 0
        self.config,result = FileLoaderDialog.getConfig()
        if result == 1:
            attrs = vars(self.config)
            print ', '.join("%s: %s" % item for item in attrs.items())
            self.analysis = Analyzer(self.config.psf,self.config.dcd, self.config.cutoff,self.config.hbondcutoff,self.config.hbondcutangle,self.config.sel1text,self.config.sel2text)
            # self.connect(self.analysis, QtCore.SIGNAL('taskUpdated'),self.handleTaskUpdated)
            # self.connect(self.analysis, QtCore.SIGNAL('frameNumberSet'),self.setFrameNumber)
            # self.progressWidget.show()
            if parallel:
                pass
            else:
                self.analysis.runFrameScan()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Data loaded: %f frames scanned." % len(self.analysis.contactResults))
            msg.setInformativeText("")
            msg.setWindowTitle("Data loaded")
            msg.setDetailedText("Now click on Analysis to proceed")
            msg.exec_()
    def handleTaskUpdated(self):
    	print self.analysis.currentFrame
    	self.progressWidget.setValue(self.analysis.currentFrame)

    def setFrameNumber(self):
    	self.progressWidget.setMax(self.analysis.totalFrameNumber)

    def analyzeParallel(self,map1,map2):
        nproc = int(self.settingsView.coreBox.value())
        start = time.time()
        trajData = self.analysis.getTrajectoryData()
        contResults = self.analysis.contactResults
        tasks = []
        results = []
        rank = 0
        manager = multiprocessing.Manager()
        d=manager.list(trajData)
        all_chunk = chunks(contResults,nproc)
        pool = multiprocessing.Pool(nproc)
        print "Running on %d cores" % nproc
        for c in all_chunk:
            results.append( pool.apply_async( loop_frame, args=(c,map1,map2,d)) )
            rank +=1
        # TODO: might be important, but without, it's faster and until now, provides the same results
        pool.close()
        pool.join()
        stop = time.time()
        print "time: ", str(stop-start), rank
        print str(len(c)), rank
        allkeys = []
        frame_contacts_accumulated = []
        print len(results)
        for res in results:
            rn = res.get()
            allkeys.extend(rn[0])
            frame_contacts_accumulated.extend(rn[1])
        accumulatedContactsDict = {}
        #   start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if not key in frame_dict:  # puts empty score TempContactAccumulate in dict
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
        # print stop - start
        glob_stop = time.time()
        print glob_stop - start
        return finalAccumulatedContacts

    def analyzeDataPushed(self):
        self.maps, result = AnalysisDialog.getMapping()
        if result == 1:
            map1 = self.maps[0]
            map2 = self.maps[1]
            nproc = int(self.settingsView.coreBox.value())
            frames = len(self.analysis.contactResults[0])
            print "frames: ",frames
            # do not allow multiprocessing unless the trajectory has enough frames
            if nproc == 1: #or  frames <= 5*nproc:
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
            # testing shit
            # maxresids1 = []
            # maxresids2 = []
            # for cont in self.contacts:
            #     cont.determineBackboneSidechainType()
            #     maxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
            #     maxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))
            # Generate some test data
            # x = np.arange(1,np.max(maxresids1)+2)
            # y = np.arange(1,np.max(maxresids2)+2)
            # data = np.zeros((len(x), len(y)))
            # for cont in self.contacts:
            #     r1 = int(cont.key1[AccumulationMapIndex.resid])
            #     r2 = int(cont.key2[AccumulationMapIndex.resid])
            #     hbonds = cont.hbondFramesScan()
            #     count = np.count_nonzero(hbonds)
            #     if count > 0:
            #         print cont.title + " contains " + str(count) + " hbonds in total"
            #     # data[r1,r2] = cont.total_time(1, 0)
            #     data[r1, r2] = cont.mean_score()
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # cax = ax.matshow(data, cmap=cm.coolwarm, aspect='equal')
            # fig.tight_layout()
            # fig.colorbar(cax)
            # plt.savefig("contactmap.png")



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
                #no weight or filters
                self.painter.showHbondScores = False
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

    def showDeveloperInfo(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        info = QLabel("Developers: Maximilian Scheurer and Peter Rodenkirch")
        info2 = QLabel("Departments: TCBG, University of Illinois at Urbana-Champaign; BZH Heidelberg University")
        mail = QLabel("Contact: mscheurer@ks.uiuc.edu, rodenkirch@stud.uni-heidelberg.de")
        copyright = QLabel("Version 0.1.1a, November 2016")

        grid.addWidget(info, 0, 0)
        grid.addWidget(info2, 1, 0)
        grid.addWidget(mail, 2, 0)
        grid.addWidget(copyright, 3, 0)

        d.setWindowTitle("Developer Info")
        d.resize(150, 80)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def pushExport(self):
        self.exportWidget = ExportTabWidget()
        self.exportWidget.valueUpdated.connect(self.handleExportUpdate)
        self.exportWidget.setContacts(self.filteredContacts)
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
        d.resize(300, 150)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def createTclScriptVis(self):
        f = open('vis.tcl', 'w')
        f.write('mol new %s \n' % self.config.psf)
        f.write('mol addfile %s \n' % self.config.dcd)
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
                currentSel1String = " and ".join(currentSel1)
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
                currentSel1String = " and ".join(currentSel1)
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



class FileLoaderDialog(QDialog):
    def __init__(self, parent = None):
        super(FileLoaderDialog, self).__init__(parent)

        self.psf = ""
        self.dcd = ""

        grid = QGridLayout(self)

        buttonPsf = QPushButton("PSF")
        buttonPsf.clicked.connect(self.pick_psf)

        buttonDcd = QPushButton("DCD")
        buttonDcd.clicked.connect(self.pick_dcd)

        grid.addWidget(buttonPsf,0,0)
        grid.addWidget(buttonDcd,0,1)

        cutoffLabel = QLabel("distance cutoff: ")
        cutoffAngleLabel = QLabel("angle cutoff: ")
        cutoffHbondLabel = QLabel("acc-h cutoff: ")
        selection1Label = QLabel("selection 1: ")
        selection2Label = QLabel("selection 2: ")

        self.cutoffField = QLineEdit("5.0")
        self.cutoffAngleField = QLineEdit("120")
        self.cutoffHbondField = QLineEdit("2.5")
        self.selection1Field = QLineEdit("segid RN11")
        self.selection2Field = QLineEdit("segid UBQ")

        grid.addWidget(cutoffLabel,1,0)
        grid.addWidget(cutoffAngleLabel,2,0)
        grid.addWidget(cutoffHbondLabel,3,0)
        grid.addWidget(selection1Label,4,0)
        grid.addWidget(selection2Label,5,0)

        grid.addWidget(self.cutoffField,1,1)
        grid.addWidget(self.cutoffAngleField,2,1)
        grid.addWidget(self.cutoffHbondField,3,1)
        grid.addWidget(self.selection1Field,4,1)
        grid.addWidget(self.selection2Field,5,1)


        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons,6,0)

    def pick_psf(self):
        psfname = QFileDialog.getOpenFileNames(self, "Open psf")
        for file in psfname[0]:
            self.psf = file
            break
    def pick_dcd(self):
        dcdname = QFileDialog.getOpenFileNames(self, "Open dcd")
        for file in dcdname[0]:
            self.dcd = file
            break

    def configuration(self):
        config = Configuration(self.psf,self.dcd, float(self.cutoffField.text()),float(self.cutoffHbondField.text()),float(self.cutoffAngleField.text()),self.selection1Field.text(),self.selection2Field.text())
        return config

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def getConfig(parent = None):
        dialog = FileLoaderDialog(parent)
        result = dialog.exec_()
        config = dialog.configuration()
        return (config, result == QDialog.Accepted)

class AnalysisDialog(QDialog):
    def __init__(self, parent = None):
        super(AnalysisDialog, self).__init__(parent)

        grid = QGridLayout(self)

        title1 = QLabel("selection 1")
        title2 = QLabel("selection 2")
        indexLabel = QLabel("index: ")
        nameLabel = QLabel("atom name: ")
        residLabel = QLabel("resid: ")
        resnameLabel = QLabel("resname: ")
        segidLabel = QLabel("segid: ")

        self.index1Checkbox = QCheckBox()
        self.name1Checkbox = QCheckBox()
        self.resid1Checkbox = QCheckBox()
        self.resname1Checkbox = QCheckBox()
        self.segid1Checkbox = QCheckBox()

        self.index2Checkbox = QCheckBox()
        self.name2Checkbox = QCheckBox()
        self.resid2Checkbox = QCheckBox()
        self.resname2Checkbox = QCheckBox()
        self.segid2Checkbox = QCheckBox()

        grid.addWidget(title1,0,1)
        grid.addWidget(title2,0,2)
        grid.addWidget(indexLabel,1,0)
        grid.addWidget(nameLabel,2,0)
        grid.addWidget(residLabel,3,0)
        grid.addWidget(resnameLabel,4,0)
        grid.addWidget(segidLabel,5,0)

        grid.addWidget(self.index1Checkbox,1,1)
        grid.addWidget(self.name1Checkbox,2,1)
        grid.addWidget(self.resid1Checkbox,3,1)
        grid.addWidget(self.resname1Checkbox,4,1)
        grid.addWidget(self.segid1Checkbox,5,1)

        grid.addWidget(self.index2Checkbox,1,2)
        grid.addWidget(self.name2Checkbox,2,2)
        grid.addWidget(self.resid2Checkbox,3,2)
        grid.addWidget(self.resname2Checkbox,4,2)
        grid.addWidget(self.segid2Checkbox,5,2)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons,6,0)


    def mapping(self):
        # atom types will not be supported in the future
        map1 = [self.index1Checkbox.isChecked(),0,self.name1Checkbox.isChecked(),self.resid1Checkbox.isChecked(),self.resname1Checkbox.isChecked(),self.segid1Checkbox.isChecked()]
        map2 = [self.index2Checkbox.isChecked(),0,self.name2Checkbox.isChecked(),self.resid2Checkbox.isChecked(),self.resname2Checkbox.isChecked(),self.segid2Checkbox.isChecked()]
        print map1, map2
        return [map1,map2]

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def getMapping(parent = None):
        dialog = AnalysisDialog(parent)
        result = dialog.exec_()
        mapping = dialog.mapping()
        return (mapping, result == QDialog.Accepted)

class Configuration(object):
    """docstring for Configuration"""
    def __init__(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        super(Configuration, self).__init__()
        self.psf = psf
        self.dcd = dcd
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.sel1text = sel1text
        self.sel2text = sel2text

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


class SettingsTabWidget(QTabWidget, Ui_settingsWindowWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QTabWidget, self).__init__(parent)
        self.setupUi(self)

class SasaWidget(QWidget, Ui_SasaWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)

        self.state = True
        self.name = None

        sip.delete(self.sasaProgressBar)
        self.sasaProgressBar = PbWidget(total=100)
        self.sasaProgressBar.setProperty("value", 0)
        self.sasaProgressBar.setTextVisible(True)
        self.sasaProgressBar.setInvertedAppearance(False)
        self.sasaProgressBar.setObjectName("sasaProgressBar")
        self.gridLayout.addWidget(self.sasaProgressBar, 4, 1, 1, 1)

        self.calcSasaButton.clicked.connect(self.calculateSasa)

        # self.totalFramesToProcess

    def calculateSasa(self):
        print "calculate SASA"
        pass

    def sasaEventListener(self):
        while self.state:
            progress = 0
            for each in sasaProgressDict.keys():
                progress += sasaProgressDict[each]
                # sasaProgressDict[each] = 0
                progress = float(progress) / float(self.totalFramesToProcess) * 100
                # update bar with delivered value
                if (101 - self.sasaProgressBar.value()) < progress:
                    self.sasaProgressBar.update_bar(101 - self.sasaProgressBar.value())
                elif progress > 0:
                    self.sasaProgressBar.update_bar(progress)

                if len(self.processes.keys()) == 0:
                    self.state = False




class PbWidget(QProgressBar):
    def __init__(self, parent=None, total=20):
        super(PbWidget, self).__init__()
        self.setMinimum(0)
        self.setMaximum(total)
        self._active = False

    def update_bar(self, to_add_number):
        while True:
            time.sleep(0.01)
            value = self.value() + to_add_number
            self.setValue(value)
            qApp.processEvents()
            if (not self._active or value >= self.maximum()):
                break
        self._active = False

    def closeEvent(self, event):
        self._active = False

class ProgessWidget(QWidget):
    def __init__(self, title):
        super(QWidget, self).__init__()
        grid = QGridLayout()
        self.title = title
        self.setGeometry(QRect(20, 40, 120, 80))
        self.progress = QProgressBar(self)
        self.progress.setGeometry(0,0,80,60)
        self.progress.setGeometry(0, 0, 80, 60)
        grid.addWidget(self.progress, 0,0)
        self.setLayout(grid)


    def setMax(self, maxi):
        self.progress.setRange(0, maxi - 1)

    def setValue(self, value):
        self.progress.setValue(value)


class ColorScheme:
    custom, bbsc = range(2)
