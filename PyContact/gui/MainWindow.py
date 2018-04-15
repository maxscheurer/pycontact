from __future__ import print_function
import multiprocessing
import warnings
import copy
import sys

import PyQt5.QtCore as QtCore
from PyQt5.QtCore import QRect, pyqtSlot, QObject, QUrl
from PyQt5.QtWidgets import (QMainWindow, QTabWidget, QLabel, QDialog, QProgressBar,
                             QApplication, QGridLayout, QFileDialog, QColorDialog, QWidget)
from PyQt5.Qt import Qt, QColor
from PyQt5.QtGui import QIntValidator
from PyQt5.QtSvg import QSvgGenerator
import numpy as np

from . import MainQtGui
from ..core.Biochemistry import vdwRadius
from .WebView import WebView
from .SasaWidgets import SasaWidget
from .MoleculeTracker import MoleculeTracker
from .Canvas import Canvas
from .Dialogues import FileLoaderDialog, AnalysisDialog
from .LoadDataDialog import LoadDataDialog
from .ExportTabWidget import ExportTabWidget
from .Statistics import Statistics
from .Plotters import *
from ..core.ContactAnalyzer_old import *
from .ErrorBox import ErrorBox
from .ErrorMessages import ErrorMessages
from ..core.LogPool import *
from ..core.aroundPatch import AroundSelection
from . import Preferences
from ..exampleData.datafiles import DEFAULTSESSION, DEFAULTSESSION_PY3
from .VMDControlPanel import VMDControlPanel
from ..core.DataHandler import DataHandler
from ..core.ContactManager import ContactManager
# from TableModels import *

multiprocessing.log_to_stderr()
np.set_printoptions(threshold=np.inf)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")


class MainWindow(QMainWindow, MainQtGui.Ui_MainWindow, QObject):
    """PyContact Application Main Window with timeline."""

    def closeEvent(self, event):
        """Closing application when Exit on MainWindow is clicked."""
        print("Closing Application")
        event.accept()
        QApplication.quit()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("PyContact")

        self.config = None
        self.maps = None
        self.contactManager = None

        # canvas contains both labels and frame boxes for drawing
        self.webView = WebView()
        self.webView.show()

        self.canvas = Canvas()
        self.scrollArea.setWidget(self.canvas)
        self.scrollArea.horizontalScrollBar().valueChanged.connect(self.horizontalScrollBarChanged)
        self.actionExportData.triggered.connect(self.pushExport)
        self.actionLoad_Data.triggered.connect(self.loadDataPushed)
        self.actionExport_Session.triggered.connect(self.exportSession)
        self.actionImport_Session.triggered.connect(self.importSession)
        self.actionShow_Info.triggered.connect(self.showDeveloperInfo)
        # settings and filters
        self.settingsView = PreferencesWidget()
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

        # apply color button, outdated?
        self.colorScheme = ColorScheme.bbsc

        self.actionDefault.triggered.connect(self.loadDefault)

        self.currentSelection1 = "-"
        self.currentSelection2 = "-"

        # setup of extra widgets
        self.exportWidget = ExportTabWidget()
        self.sasaView = SasaWidget()
        self.statisticsView = None

        self.moleculeTracker = MoleculeTracker()
        self.actionTrack_Molecule.triggered.connect(self.showMoleculeTracker)

        self.vismode = False
        self.visModeButton.setCheckable(True)
        self.visModeButton.setChecked(False)
        self.visModeButton.clicked.connect(self.switchedToVisMode)

        self.vmdpanel = VMDControlPanel()
        self.actionVMD_Remote_Control.triggered.connect(self.showVMDControlPanel)

        self.canvas.clickedRowSignal.connect(self.updateVMDSelections)
        self.canvas.clickedColumnSignal.connect(self.updateVMDFrame)
        self.updateSettings()
        # TODO: enable
        # self.updateFilters()

        self.actionDefault.setText("Load sample data")

    def horizontalScrollBarChanged(self):
        x = self.scrollArea.horizontalScrollBar().value()
        y = self.canvas.labelView.y()
        self.canvas.labelView.move(x, y)

    def showVMDControlPanel(self):
        """Shows the VMD control panel, to remotely access VMD from PyContact."""
        self.vmdpanel.show()

    def showMoleculeTracker(self):
        """Shows the VMD control panel, to remotely access VMD from PyContact."""
        self.moleculeTracker.show()

    def showContactAreaView(self):
        """Shows the SASA computation panel."""
        self.sasaView.nsPerFrame = float(self.settingsView.nsPerFrameField.text())
        self.sasaView.show()
        # TODO: ask user which trajectory should be used for SASA
        # if self.contactManager:
        #     self.sasaView.setFilePaths(self.contactManager.getFilePaths())

    def switchedToVisMode(self):
        """Switch to vis mode, to show selected contacts directly in VMD."""
        if self.visModeButton.isChecked():
            self.vismode = True
            # conversions with clicked frames are not allowed
            self.frameStrideField.setText("1")
        else:
            self.vismode = False
        self.canvas.switchToVisMode(self.vismode)
        self.updateSettings()
        self.updateFilters()

    @pyqtSlot()
    def updateVMDSelections(self):
        """Updates the selected contact in VMD via the vmd panel."""
        # TODO: refactor for contactManager
        if self.vmdpanel.connected:
            self.vmdpanel.updateSelections(
                self.contactManager.trajectoryScanParameters.sel1text,
                self.contactManager.trajectoryScanParameters.sel2text,
                [])

    @pyqtSlot()
    def updateVMDFrame(self):
        """Updates the selected frame in VMD via the vmd panel."""
        if self.vmdpanel.connected:
            self.vmdpanel.gotoVMDFrame(self.canvas.clickedColumn)

    def updateSelectionLabels(self, sel1, sel2):
        """Updates the current selection in the info labels."""
        self.currentSelection1 = sel1
        self.currentSelection2 = sel2
        self.selection1label.setText(sel1)
        self.selection2label.setText(sel2)

    def importSession(self):
        """Imports a saved session from file."""
        fnames = QFileDialog.getOpenFileNames(self, "Open file")
        importfile = ""
        for f in fnames[0]:
            importfile = f
            break
        if importfile == "" or len(fnames) == 0:
            return

        # TODO: refactor with contactManager
        # self.contacts, arguments, trajArgs, self.maps, contactResults = DataHandler.importSessionFromFile(importfile)
        # self.analysis = Analyzer(*arguments)
        # self.analysis.contactResults = contactResults
        # self.analysis.setTrajectoryData(*trajArgs)
        # self.analysis.finalAccumulatedContacts = self.contacts
        # self.sasaView.setFilePaths(*self.analysis.getFilePaths())
        # self.exportWidget.setFilePaths(*self.analysis.getFilePaths())
        # self.updateSelectionLabels(arguments[5], arguments[6])
        # self.updateSettings()
        # self.updateFilters()

    def exportSession(self):
        """Exports the current session to file."""
        fileName = QFileDialog.getSaveFileName(self, 'Export file')
        filestring = fileName[0]
        if filestring == "":
            return
        # TODO: refactor with contactManager
        # if self.contacts is not None and self.analysis is not None:
        #     self.setInfoLabel("Exporting current session...")
        #     DataHandler.writeSessionToFile(filestring, self.analysis)
        #     self.cleanInfoLabel()
        # else:
        #     box = ErrorBox(ErrorMessages.NOEXPDATA)
        #     box.exec_()
        #     return

    def loadDefault(self):
        """Loads the default session."""

        # TODO: refactor with contactManager
        # if (sys.version_info > (3, 0)):
        #     self.contacts, arguments, trajArgs, self.maps, contactResults = DataHandler.importSessionFromFile(DEFAULTSESSION_PY3)
        # else:
        #     self.contacts, arguments, trajArgs, self.maps, contactResults = DataHandler.importSessionFromFile(DEFAULTSESSION)
        # self.analysis = Analyzer(*arguments)
        # self.analysis.contactResults = contactResults
        # self.analysis.setTrajectoryData(*trajArgs)
        # self.analysis.finalAccumulatedContacts = self.contacts
        # self.sasaView.setFilePaths(*self.analysis.getFilePaths())
        # self.exportWidget.setFilePaths(*self.analysis.getFilePaths())
        # self.updateSelectionLabels(arguments[5], arguments[6])
        # self.updateSettings()
        # self.updateFilters()

    def loadDataPushed(self):
        """Loads the trajectory data with the chosen initial parameters."""
        self.config, result = LoadDataDialog.getConfig()
        if result == 1:
            QApplication.processEvents()
            self.setInfoLabel("Loading trajectory and running atomic contact analysis...")
            nproc = int(self.settingsView.coreBox.value())
            #self.analysis = Analyzer(self.config.psf, self.config.dcd, self.config.cutoff, self.config.hbondcutoff,
                                     # self.config.hbondcutangle, self.config.sel1text, self.config.sel2text)
            self.contactManager = ContactManager(self.config.topology,
                                                 self.config.trajectories,
                                                 self.config.cutoff,
                                                 self.config.hbondcutoff,
                                                 self.config.hbondcutangle,
                                                 self.config.sel1text,
                                                 self.config.sel2text)
            QApplication.processEvents()
            # TODO: enable
            # try:
                # self.analysis.runFrameScan(nproc)
            self.contactManager.readTrajectories(nproc)
            # except:
            #     box = ErrorBox("Error while loading data: Probably you specified an atom selection with 0 atoms or invalid input files.")
            #     box.exec_()
            #     self.loadDataPushed()

            # self.setInfoLabel("%d frames loaded." % len(self.analysis.contactResults))
            self.updateSelectionLabels(self.config.sel1text, self.config.sel2text)
            # self.sasaView.setFilePaths(*self.analysis.getFilePaths())
            # self.exportWidget.setFilePaths(*self.analysis.getFilePaths())

    @pyqtSlot(float)
    def updateAnalyzedFrames(self, value):
        """Handles the progress bar update."""
        # print("Updating frames", value)
        self.progressBar.setValue(100 * value)
        QApplication.processEvents()

    def setInfoLabel(self, txt):
        """Sets the Info label text."""
        self.statusLabel.setText(txt)

    def cleanInfoLabel(self):
        """Clears the Info label text."""
        self.setInfoLabel("-")

    def analyzeDataPushed(self):
        # TODO: refactor
        """Handles the Analyzer after the Accumulation maps have been set."""
        # if self.analysis is None:
        #     box = ErrorBox(ErrorMessages.NODATA_PROMPTLOAD)
        #     box.exec_()
        #     return

        self.maps, result = AnalysisDialog.getMapping()
        if result == 1:
            # self.analysis.frameUpdate.connect(self.updateAnalyzedFrames)
            self.setInfoLabel("Analyzing contacts...")
            self.contactManager.accumulateContacts(*self.maps)
            self.canvas.setAccumulatedTrajectory(self.contactManager.accumulatedContactTrajectories[0])
            self.updateCanvas()
            # self.contacts = self.analysis.runContactAnalysis(map1, map2, nproc)
            # self.progressBar.setValue(0)
            # self.setInfoLabel("Updating timeline...")
            QApplication.processEvents()
            # self.updateSettings()
            # self.updateFilters()
            # self.cleanInfoLabel()

    def updateCanvas(self):
        print("update Canvas")
        self.canvas.rendered = False
        self.canvas.repaint()
        self.canvas.update()

    def updateSettings(self):
        """Updates the settings chosen from the settings view."""
        self.canvas.nsPerFrame = float(self.settingsView.nsPerFrameField.text())
        self.canvas.threshold = float(self.settingsView.thresholdField.text())
        self.canvas.rendered = False
        self.canvas.colorScheme = self.colorScheme
        self.canvas.customColor = self.customColor
        self.canvas.repaint()
        self.canvas.update()
        self.sasaView.nsPerFrame = float(self.settingsView.nsPerFrameField.text())

    def updateFilters(self):
        """Updates the chosen filters in MainWindow."""
        if self.vismode is True:
            self.frameStrideField.setText("1")
        stride = int(self.frameStrideField.text())
        if stride < 1:
            stride = 1
            QApplication.processEvents()
            self.frameStrideField.setText(str(stride))
        # print("stride: ", stride)
        self.canvas.merge = stride
        self.canvas.labelView.clean()
        self.canvas.showHbondScores = False
        # total time filter
        totalTimeActive = self.activeTotalTimeCheckbox.isChecked()
        scoreActive = self.activeScoreCheckbox.isChecked()
        sortingActive = self.activeSortingBox.isChecked()
        onlyActive = self.onlyBoxActiveCheckbox.isChecked()
        filterActive = (totalTimeActive or scoreActive or sortingActive or onlyActive)
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
            self.canvas.range = [lower, upper]
            self.canvas.rangeFilterActive = False
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
                self.canvas.rangeFilterActive = True
                frameRangeFilter = FrameFilter("framer")
                self.filteredContacts = frameRangeFilter.extractFrameRange(self.filteredContacts, [lower, upper])
            for c in self.filteredContacts:
                c.setScores()
                c.setContactType()
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
                        sorter.setThresholdAndNsPerFrame(float(self.settingsView.thresholdField.text()),
                                                         float(self.settingsView.nsPerFrameField.text()))
                        self.filteredContacts = sorter.sortContacts(self.filteredContacts)
                    if onlyActive:
                        key = self.selectOnlyToolbox.currentText()
                        only = OnlyFilter("only", key, 0)
                        self.filteredContacts = only.filterContacts(self.filteredContacts)
                        if key == "hbonds":
                            self.canvas.showHbondScores = True
                    self.canvas.contacts = self.filteredContacts
                    self.canvas.rendered = False
                    self.canvas.repaint()
                    self.canvas.update()
                    if len(self.filteredContacts) == 0:
                        self.canvas.labelView.clean()
            else:
                # no weight or filters
                self.canvas.showHbondScores = False
                self.canvas.contacts = self.filteredContacts
                self.canvas.rendered = False
                self.canvas.repaint()
                self.canvas.update()

        # Update data for export
        self.exportWidget.setContacts(self.filteredContacts)
        if self.maps is not None:
            self.exportWidget.setMaps(self.maps[0], self.maps[1])
            self.exportWidget.setMapLabels(self.analysis.sel1text, self.analysis.sel2text)
            self.vmdpanel.sel1 = self.analysis.sel1text
            self.vmdpanel.sel2 = self.analysis.sel2text
            self.vmdpanel.filteredContactList = self.filteredContacts
        self.exportWidget.setThresholdAndNsPerFrame(self.canvas.threshold, self.canvas.nsPerFrame)

    def openPrefs(self):
        """Opens the preferences panel."""
        self.settingsView.show()

    def showStatistics(self):
        """Shows general statistics of the analyzed data over all frames."""
        # TODO: refactor with contactManager
        # if len(self.contacts) == 0 or self.contacts is None:
        #     box = ErrorBox(ErrorMessages.NOSCORES_PROMPTANALYSIS)
        #     box.exec_()
        #     return
        # self.statisticsView = Statistics(self.contacts, float(self.settingsView.nsPerFrameField.text()))
        # self.statisticsView.showNormal()

    def showDeveloperInfo(self):
        """Shows information about the contributing authors."""
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        info = QLabel("Developers: Maximilian Scheurer and Peter Rodenkirch")
        info2 = QLabel("Departments: TCBG, University of Illinois at Urbana-Champaign; BZH Heidelberg University")
        mail = QLabel("Contact: mscheurer@ks.uiuc.edu, rodenkirch@stud.uni-heidelberg.de")
        copyright = QLabel("Version 1.0, May 2017")

        grid.addWidget(info, 0, 0)
        grid.addWidget(info2, 1, 0)
        grid.addWidget(mail, 2, 0)
        grid.addWidget(copyright, 3, 0)

        d.setWindowTitle("Developer Info")
        d.resize(150, 80)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def pushExport(self):
        """Opens the export panel."""
        # TODO: refactor with contactManager
        # self.exportWidget.valueUpdated.connect(self.handleExportUpdate)
        # self.exportWidget.setContacts(self.filteredContacts)
        # if self.maps is not None:
        #     self.exportWidget.setMaps(self.maps[0], self.maps[1])
        #     self.exportWidget.setMapLabels(self.analysis.sel1text, self.analysis.sel2text)
        # self.exportWidget.setThresholdAndNsPerFrame(self.canvas.threshold, self.canvas.nsPerFrame)
        # self.exportWidget.show()

    @QtCore.Slot(str, str)
    def handleExportUpdate(self, fileName, fileType):
        """Handles the paint event after the export of the current view has been initiated."""
        print("test")
        if fileType == "PNG":
            if len(fileName) > 0:
                print("Saving current view to ", fileName)
                currentView = self.canvas.grab()
                currentView.save(fileName)
        elif fileType == "SVG":
            if len(fileName) > 0:
                print("Saving current view to ", fileName)
                generator = QSvgGenerator()
                generator.setFileName(fileName)
                generator.setSize(self.canvas.size())
                generator.setViewBox(self.canvas.rect())
                self.canvas.renderContact(generator)
        self.canvas.rendered = False
        self.canvas.repaint()
        self.canvas.update()

    def showColorPicker(self):
        """Shows a color picker for the current view."""
        col = QColorDialog.getColor()
        self.customColor = col
        if col.isValid():
            self.settingsView.pickColorButton.setStyleSheet("QWidget { background-color: %s }" %
                                                            self.customColor.name())


class PreferencesWidget(QTabWidget, Preferences.Ui_PreferencesPanel):
    """Defines the preferences panel"""
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)
        self.setupUi(self)


class ColorScheme:
    custom, bbsc = range(2)
