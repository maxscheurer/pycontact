from __future__ import print_function

from PyQt5.QtWidgets import QDialog, QGridLayout, QPushButton, QLabel, QLineEdit, QDialogButtonBox, QFileDialog, QCheckBox
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator

from ..core.LoadConfiguration import Configuration
from PyContact.core.ContactAnalyzer import Analyzer
from HelpButton import HelpButton


class TopoTrajLoaderDialog(QDialog):
    """Dialog to load the topology and trajectory file."""
    def __init__(self, parent=None):
        super(TopoTrajLoaderDialog, self).__init__(parent)
        self.setWindowTitle("Load Data")
        self.psf = ""
        self.dcd = ""

        grid = QGridLayout(self)

        buttonPsf = QPushButton("Topology")
        buttonPsf.clicked.connect(self.pick_psf)

        buttonDcd = QPushButton("Trajectory")
        buttonDcd.clicked.connect(self.pick_dcd)

        helpButton = HelpButton()

        grid.addWidget(buttonPsf, 0, 0)
        grid.addWidget(buttonDcd, 0, 1)
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons, 1, 0)
        grid.addWidget(helpButton, 2, 0)

    def pick_psf(self):
        """Pick the topology file."""
        psfname = QFileDialog.getOpenFileNames(self, "Open topology")
        for file in psfname[0]:
            self.psf = file
            break

    def pick_dcd(self):
        """Pick the trajectory file."""
        dcdname = QFileDialog.getOpenFileNames(self, "Open trajectory")
        for file in dcdname[0]:
            self.dcd = file
            break

    def configuration(self):
        """Returns the chosen configuration."""
        config = [self.psf, self.dcd]
        return config

    @staticmethod
    def getConfig(parent=None):
        """Static method to create the dialog and return (date, time, accepted)."""
        dialog = TopoTrajLoaderDialog(parent)
        result = dialog.exec_()
        config = dialog.configuration()
        return config, result == QDialog.Accepted


class FileLoaderDialog(QDialog):
    """Initial file loader dialog, including initial parameter settings."""
    def __init__(self, parent=None):
        super(FileLoaderDialog, self).__init__(parent)
        self.setWindowTitle("Load Data")
        self.psf = ""
        self.dcd = ""
        production = 1

        grid = QGridLayout(self)

        buttonPsf = QPushButton("Topology")
        buttonPsf.clicked.connect(self.pick_psf)

        buttonDcd = QPushButton("Trajectory")
        buttonDcd.clicked.connect(self.pick_dcd)

        grid.addWidget(buttonPsf,0,0)
        grid.addWidget(buttonDcd,0,1)

        cutoffLabel = QLabel("distance cutoff: ")
        cutoffAngleLabel = QLabel("angle cutoff: ")
        cutoffHbondLabel = QLabel("acc-h cutoff: ")
        selection1Label = QLabel("selection 1: ")
        selection2Label = QLabel("selection 2: ")

        self.cutoffField = QLineEdit("5.0")
        posDoubleValidator = QDoubleValidator()
        posDoubleValidator.setBottom(0)
        self.cutoffField.setValidator(posDoubleValidator)
        self.cutoffAngleField = QLineEdit("120")
        self.cutoffAngleField.setValidator(posDoubleValidator)
        self.cutoffHbondField = QLineEdit("2.5")
        self.cutoffHbondField.setValidator(posDoubleValidator)

        if production:
            self.selection1Field = QLineEdit("")
            self.selection2Field = QLineEdit("")
        else:
            self.selection1Field = QLineEdit("segid RN11")
            self.selection2Field = QLineEdit("segid UBQ")

        grid.addWidget(cutoffLabel, 1, 0)
        grid.addWidget(cutoffAngleLabel, 2, 0)
        grid.addWidget(cutoffHbondLabel, 3, 0)
        grid.addWidget(selection1Label, 4, 0)
        grid.addWidget(selection2Label, 5, 0)

        grid.addWidget(self.cutoffField, 1, 1)
        grid.addWidget(self.cutoffAngleField, 2, 1)
        grid.addWidget(self.cutoffHbondField, 3, 1)
        grid.addWidget(self.selection1Field, 4, 1)
        grid.addWidget(self.selection2Field, 5, 1)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons, 6, 0)

    def pick_psf(self):
        """Pick topology file."""
        psfname = QFileDialog.getOpenFileNames(self, "Open topology")
        for file in psfname[0]:
            self.psf = file
            break

    def pick_dcd(self):
        """Pick trajectory file."""
        dcdname = QFileDialog.getOpenFileNames(self, "Open trajectory")
        for file in dcdname[0]:
            self.dcd = file
            break

    def configuration(self):
        """Returns the chosen configuration."""
        config = Configuration(self.psf, self.dcd, float(self.cutoffField.text()), float(self.cutoffHbondField.text()),
                               float(self.cutoffAngleField.text()), self.selection1Field.text(),
                               self.selection2Field.text())
        return config

    @staticmethod
    def getConfig(parent=None):
        """Static method to create the dialog and return (date, time, accepted)."""
        dialog = FileLoaderDialog(parent)
        result = dialog.exec_()
        config = dialog.configuration()
        return config, result == QDialog.Accepted


class AnalysisDialog(QDialog):
    """Dialog to define the Accumulation maps."""
    def __init__(self, parent=None):
        super(AnalysisDialog, self).__init__(parent)

        grid = QGridLayout(self)

        self.title1 = QLabel("selection 1")
        self.title2 = QLabel("selection 2")
        indexLabel = QLabel("index: ")
        nameLabel = QLabel("atom name: ")
        residLabel = QLabel("resid: ")
        resnameLabel = QLabel("resname: ")
        segidLabel = QLabel("segid: ")

        self.setWindowTitle("Analysis - Score Accumulation")

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

        grid.addWidget(self.title1, 0, 1)
        grid.addWidget(self.title2, 0, 2)
        grid.addWidget(indexLabel, 1, 0)
        grid.addWidget(nameLabel, 2, 0)
        grid.addWidget(residLabel, 3, 0)
        grid.addWidget(resnameLabel, 4, 0)
        grid.addWidget(segidLabel, 5, 0)

        grid.addWidget(self.index1Checkbox, 1, 1)
        grid.addWidget(self.name1Checkbox, 2, 1)
        grid.addWidget(self.resid1Checkbox, 3, 1)
        grid.addWidget(self.resname1Checkbox, 4, 1)
        grid.addWidget(self.segid1Checkbox, 5, 1)

        grid.addWidget(self.index2Checkbox, 1, 2)
        grid.addWidget(self.name2Checkbox, 2, 2)
        grid.addWidget(self.resid2Checkbox, 3, 2)
        grid.addWidget(self.resname2Checkbox, 4, 2)
        grid.addWidget(self.segid2Checkbox, 5, 2)

        # OK and Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        grid.addWidget(buttons, 6, 0)
        self.gridLayout = grid

    def mapping(self):
        """Creates the Accumulation maps from the checkbox values."""
        # atom types will not be supported in the future
        map1 = [self.index1Checkbox.isChecked(), self.name1Checkbox.isChecked(), self.resid1Checkbox.isChecked(),
                self.resname1Checkbox.isChecked(), self.segid1Checkbox.isChecked()]
        map2 = [self.index2Checkbox.isChecked(), self.name2Checkbox.isChecked(), self.resid2Checkbox.isChecked(),
                self.resname2Checkbox.isChecked(), self.segid2Checkbox.isChecked()]
        print("Accumulation maps: ", map1, map2)
        return [map1, map2]

    @staticmethod
    def getMapping(parent=None):
        """Static method to create the dialog and return (date, time, accepted)."""
        dialog = AnalysisDialog(parent)
        result = dialog.exec_()
        mapping = dialog.mapping()
        return mapping, result == QDialog.Accepted


class AnalysisSingleDialog(AnalysisDialog):
    def __init__(self, parent=None):
        super(AnalysisSingleDialog, self).__init__(parent)
        self.index2Checkbox.setHidden(True)
        self.name2Checkbox.setHidden(True)
        self.resid2Checkbox.setHidden(True)
        self.resname2Checkbox.setHidden(True)
        self.segid2Checkbox.setHidden(True)
        self.title2.setHidden(True)
        self.title1.setText("selection")
        self.setWindowTitle("Molecule Tracking Selection")

    @staticmethod
    def getMapping(parent=None):
        """Static method to create the dialog and return (date, time, accepted)."""
        dialog = AnalysisSingleDialog(parent)
        result = dialog.exec_()
        mapping = dialog.mapping()
        return mapping[0], result == QDialog.Accepted
