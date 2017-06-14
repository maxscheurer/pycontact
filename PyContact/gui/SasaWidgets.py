from __future__ import print_function
import sip
import time
import os


from PyQt5.QtWidgets import QWidget, QProgressBar, QApplication, QFileDialog
import MDAnalysis
import numpy as np

from .Plotters import SimplePlotter
from sasa_gui import *
from ..core.multi_accumulation import chunks
from ..core.Biochemistry import vdwRadius
from ..core.LogPool import *
from ..cy_modules import cy_gridsearch
from Dialogues import TopoTrajLoaderDialog
from ErrorBox import ErrorBox
from ErrorMessages import ErrorMessages

# manage processes for SASA
sasaProgressManager = multiprocessing.Manager()
sasaProgressDict = sasaProgressManager.dict()
np.set_printoptions(threshold=np.inf)


def calculate_sasa_parallel(input_coords, natoms, pairdist, nprad,
                            surfacePoints, probeRadius, pointstyle,
                            restricted, restrictedList, rank):
    """Computes the SASA in parallel."""

    temp_sasa = []
    frames_processed = 0
    sasaProgressDict[rank] = frames_processed
    # print(len(input_coords))
    for c in input_coords:
        coords = np.reshape(c, (1, natoms * 3))
        npcoords = np.array(coords, dtype=np.float32)
        # print("start C")
        # startC = time.time()
        asa = cy_gridsearch.cy_sasa(npcoords, natoms, pairdist, 0, -1, nprad, surfacePoints, probeRadius,
                                    pointstyle, restricted, restrictedList)
        # stopC = time.time()
        # print("time for grid search: ", (stopC - startC))
        # print("asa:", asa)
        temp_sasa.append(asa)
        frames_processed += 1
        sasaProgressDict[rank] = frames_processed
    return temp_sasa


class SasaWidget(QWidget, Ui_SasaWidget):
    """Provides a UI for the SASA calculation, including restriction selection and contact area calculation."""
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)

        self.state = True
        self.name = None
        self.psf, self. dcd = "", ""
        self.allSasas = []
        self.totalFramesToProcess = 0
        self.nsPerFrame = 1.0

        sip.delete(self.sasaProgressBar)
        self.sasaProgressBar = PbWidget(total=100)
        self.sasaProgressBar.setProperty("value", 0)
        self.sasaProgressBar.setTextVisible(True)
        self.sasaProgressBar.setInvertedAppearance(False)
        self.sasaProgressBar.setObjectName("sasaProgressBar")
        self.previewPlot = SimplePlotter(None, width=4, height=2, dpi=70)
        self.graphGridLayout.addWidget(self.previewPlot)
        self.gridLayout.addWidget(self.sasaProgressBar, 8, 1, 1, 2)
        self.calcSasaButton.clicked.connect(self.calculateSasa)
        self.loadDataButton.clicked.connect(self.loadData)
        self.clearDataButton.clicked.connect(self.clearData)
        self.savePlotButton.clicked.connect(self.savePlot)
        self.exportDataButton.clicked.connect(self.exportData)
        self.topoloader = TopoTrajLoaderDialog()

    def setFilePaths(self, *argv):
        """Sets the current trajectory paths from the main view."""
        self.psf = argv[0][0]
        self.dcd = argv[0][1]

    def loadData(self):
        """Sets the chosen trajectory and topology paths."""
        loadedData = self.topoloader.getConfig()
        self.psf = loadedData[0][0]
        self.dcd = loadedData[0][1]

    def clearData(self):
        """Clears the whole view and sets it to the initial state."""
        self.psf = ""
        self.dcd = ""
        self.allSasas = []
        self.totalFramesToProcess = 0
        self.sasaProgressBar.setProperty("value", 0)
        sip.delete(self.previewPlot)
        self.previewPlot = SimplePlotter(None, width=4, height=2, dpi=70)
        self.graphGridLayout.addWidget(self.previewPlot)
        self.sasaSelection1TextField.setText("")
        self.sasaSelection2TextField.setText("")
        self.sasaRestrictionTextField.setText("")


    def savePlot(self):
        """Saves the generated SASA plot over all frames."""
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            if file_extension == "":
                file_extension = ".png"
            file_extension = file_extension[1:]
            try:
                self.previewPlot.saveFigure(path, file_extension)
            except ValueError:
                box = ErrorBox("File format " + file_extension + " is not supported.\nPlease choose from eps, pdf, pgf,"
                                                                 " png, ps, raw, rgba, svg, svgz. ")
                box.exec_()

    def exportData(self):
        """Exports the computed SASA values of each frame to a text file."""
        fileName = QFileDialog.getSaveFileName(self, 'Export Path')
        if len(fileName[0]) > 0:
            path, file_extension = os.path.splitext(fileName[0])
            if file_extension == "":
                file_extension = ".dat"

            f = open(path + file_extension, "w")
            for i in range(self.totalFramesToProcess):
                f.write(str(i) + "\t" + str(self.allSasas[i]) + "\n")
            f.close()

    def calculateSasa(self):
        """Computes the SASA of the given selections."""
        print("calculate SASA")

        self.allSasas = []

        # load psf and trajectory, make lists with radii and coordinates
        if self.psf == "" or self.dcd == "":
            e = ErrorBox(ErrorMessages.CHOOSEFILE)
            e.exec_()
            return

        try:
            u = MDAnalysis.Universe(self.psf, self.dcd)
        except IOError:
            e = ErrorBox(ErrorMessages.FILE_NOT_FOUND)
            e.exec_()
            return

        probeRadius = 1.4

        # seltext = "segid UBQ"
        # resseltext = "segid UBQ and same residue as around 5.0 (segid RN11)"
        seltext = self.sasaSelection1TextField.text()
        seltext2 = self.sasaSelection2TextField.text()
        resseltext = self.sasaRestrictionTextField.text()

        # 0=spiral, 1=random (VMD)
        pointstyle = 1
        # number of points to approximate the sphere
        surfacePoints = 50
        # pair distance
        pairdist = 2.0 * (2.0 + 1.4)

        if resseltext != "":
            restricted = 1
        else:
            restricted = 0

        selection = u.select_atoms(seltext)

        # natoms = len(selection.atoms)
        radius = []
        restrictedList = []
        if restricted:
            ressel = u.select_atoms(resseltext)
            for s in selection.atoms:
                if s in ressel.atoms:
                    restrictedList.append(1)
                else:
                    restrictedList.append(0)
                radius.append(vdwRadius(s.name[0]))
        else:
            restrictedList = [0]
            for s in selection.atoms:
                radius.append(vdwRadius(s.name[0]))
        natoms = len(selection)
        nprad = np.array(radius, dtype=np.float32)
        restrictedList = np.array(restrictedList, dtype=np.int32)

        # TODO: bug if selection is not static for all frames
        # TODO: dynamic allocation of positions in every frame!
        input_coords = []
        for ts in u.trajectory:
            # ressel = u.select_atoms(resseltext)
            # print("restricted: ", len(ressel.atoms))
            input_coords.append(selection.positions)

        nprocs = self.coreBox.value()
        input_chunks = chunks(input_coords, nprocs)
        pool = LoggingPool(nprocs)
        results = []
        rank = 0
        trajLength = len(u.trajectory)
        self.totalFramesToProcess = trajLength
        for input_coords_chunk in input_chunks:
            results.append(pool.apply_async(calculate_sasa_parallel, args=(input_coords_chunk, natoms, pairdist, nprad,
                                                                           surfacePoints, probeRadius, pointstyle,
                                                                           restricted, restrictedList, rank)))
            rank += 1
        print("ranks", rank)
        self.state = True
        self.sasaEventListener()
        pool.close()
        pool.join()

        self.state = False
        for r in results:
            self.allSasas.extend(r.get())

        del radius

        if self.calculateContactAreaCheckbox.isChecked():
            print("Calculate contact area")
            selection2 = u.select_atoms(seltext2)

            # natoms2 = len(selection2.atoms)
            radius2 = []
            restrictedList2 = []
            if restricted:
                ressel = u.select_atoms(resseltext)
                for s in selection2.atoms:
                    if s in ressel.atoms:
                        restrictedList2.append(1)
                    else:
                        restrictedList2.append(0)
                    radius2.append(vdwRadius(s.name[0]))
            else:
                print("You need a restricted selection for contact areas!")

            natoms2 = len(selection2)
            nprad = np.array(radius2, dtype=np.float32)
            restrictedList2 = np.array(restrictedList2, dtype=np.int32)

            input_coords2 = []
            for ts in u.trajectory:
                input_coords2.append(selection2.positions)

            input_chunks2 = chunks(input_coords2, nprocs)
            pool = LoggingPool(nprocs)
            results = []
            rank = 0
            trajLength = len(u.trajectory)
            self.totalFramesToProcess = trajLength
            for input_coords_chunk2 in input_chunks2:
                results.append(pool.apply_async(calculate_sasa_parallel, args=(input_coords_chunk2, natoms2, pairdist,
                                                                               nprad, surfacePoints, probeRadius,
                                                                               pointstyle, restricted, restrictedList2,
                                                                               rank)))
                rank += 1
            print("ranks", rank)
            self.state = True
            self.sasaEventListener()
            pool.close()
            pool.join()

            all_sasas2 = []
            for r in results:
                all_sasas2.extend(r.get())

            diff_list = []
            for sasa_value1, sasa_value2 in zip(self.allSasas, all_sasas2):
                diff_list.append(sasa_value1 - sasa_value2)
            self.allSasas = diff_list

        sip.delete(self.previewPlot)
        self.previewPlot = SimplePlotter(None, width=4, height=2, dpi=70)
        self.previewPlot.plot(np.arange(0, trajLength, 1) * self.nsPerFrame, self.allSasas)
        self.previewPlot.axes.set_xlabel("time [ns]")
        if self.calculateContactAreaCheckbox.isChecked():
            self.previewPlot.axes.set_ylabel(r'Contact Area [A$^{\circ}$$^{2}$]')
        else:
            self.previewPlot.axes.set_ylabel(r'SASA [A$^{\circ}$$^{2}$]')
        self.graphGridLayout.addWidget(self.previewPlot)
        self.previewPlot.update()

    def sasaEventListener(self):
        """Event listener for progress bar updates."""
        while self.state:
            QApplication.processEvents()
            progress = 0
            for each in sasaProgressDict.keys():
                progress += sasaProgressDict[each]
                # sasaProgressDict[each] = 0
            progress = float(progress) / float(self.totalFramesToProcess) * 100
            # if (101 - self.sasaProgressBar.value()) < progress:
            #     self.sasaProgressBar.update_bar(101 - self.sasaProgressBar.value())
            if progress > 0:
                self.sasaProgressBar.setValue(progress)

            if int(progress) == 100:
                # print("finished")
                for each in sasaProgressDict.keys():
                    sasaProgressDict[each] = 0
                # progress = 0
                self.state = False


class PbWidget(QProgressBar):
    """Subclassed progressbar for the SASA UI."""
    def __init__(self, total=100):
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
            if not self._active or value >= self.maximum():
                break
        self._active = False

    def closeEvent(self, event):
        self._active = False
