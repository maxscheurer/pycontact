from __future__ import print_function
import sip
import time

from PyQt5.QtWidgets import QWidget, QProgressBar
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
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.setupUi(self)

        self.state = True
        self.name = None
        self.psf, self. dcd = "", ""

        sip.delete(self.sasaProgressBar)
        self.sasaProgressBar = PbWidget(total=100)
        self.sasaProgressBar.setProperty("value", 0)
        self.sasaProgressBar.setTextVisible(True)
        self.sasaProgressBar.setInvertedAppearance(False)
        self.sasaProgressBar.setObjectName("sasaProgressBar")
        self.previewPlot = SimplePlotter(None, width=5, height=2, dpi=60)
        self.graphGridLayout.addWidget(self.previewPlot)
        self.gridLayout.addWidget(self.sasaProgressBar, 7, 1, 1, 1)
        self.calcSasaButton.clicked.connect(self.calculateSasa)
        self.loadDataButton.clicked.connect(self.loadData)
        self.topoloader = TopoTrajLoaderDialog()

    def loadData(self):
        loadedData = self.topoloader.getConfig()
        self.psf = loadedData[0][0]
        self.dcd = loadedData[0][1]

    def calculateSasa(self):
        print("calculate SASA")

        # test data:
        # psf = "rpn11_ubq_interface-ionized.psf"
        # pdb = "rpn11_ubq_interface-ionized.pdb"
        # dcd = "/home/max/Projects/pycontact/short.dcd"

        # Rafaels COH3 DOC3
        # psf = "/mnt/workspace/pycontactData/nowater.psf"
        # dcd = "/mnt/workspace/pycontactData/trajectory_short.dcd"

        # load psf and trajectory, make lists with radii and coordinates
        if self.psf == "" or self.dcd == "":
            e = ErrorBox(ErrorMessages.CHOOSEFILE)
            e.exec_()
            return

        u = MDAnalysis.Universe(self.psf, self.dcd)

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

        # not used at the moment
        # if perres:
        #     resids = sorted(set(selection.resids))
        #     segs = sorted(set(selection.segids))
        # else:
        #     pass

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
        all_sasas = []
        for r in results:
            all_sasas.extend(r.get())

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
            for sasa_value1, sasa_value2 in zip(all_sasas, all_sasas2):
                diff_list.append(sasa_value1 - sasa_value2)
            all_sasas = diff_list

        sip.delete(self.previewPlot)
        self.previewPlot = SimplePlotter(None, width=5, height=2, dpi=60)
        self.previewPlot.plot(np.arange(0, trajLength, 1), all_sasas)
        self.previewPlot.axes.set_xlabel("frame")
        self.previewPlot.axes.set_ylabel(r'SASA [A$^\circ$^2]')
        self.graphGridLayout.addWidget(self.previewPlot)
        self.previewPlot.update()

    def sasaEventListener(self):
        while self.state:
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
