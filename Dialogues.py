from __future__ import print_function
from PyQt5.QtWidgets import QDialog, QGridLayout, QPushButton, QLabel,QLineEdit,QDialogButtonBox, QFileDialog, QCheckBox
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QDoubleValidator, QRegExpValidator
from LoadConfiguration import Configuration
from PyQt5.QtCore import QRegExp
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
        posDoubleValidator = QDoubleValidator()
        posDoubleValidator.setBottom(0)
        self.cutoffField.setValidator(posDoubleValidator)
        self.cutoffAngleField = QLineEdit("120")
        self.cutoffAngleField.setValidator(posDoubleValidator)
        self.cutoffHbondField = QLineEdit("2.5")
        self.cutoffHbondField.setValidator(posDoubleValidator)

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
        print(map1, map2)
        return [map1,map2]

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def getMapping(parent = None):
        dialog = AnalysisDialog(parent)
        result = dialog.exec_()
        mapping = dialog.mapping()
        return (mapping, result == QDialog.Accepted)
