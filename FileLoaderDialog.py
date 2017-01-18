from PyQt5.QtWidgets import QDialog, QGridLayout, QPushButton, QLabel,QLineEdit,QDialogButtonBox, QFileDialog
from PyQt5.QtCore import Qt

from LoadConfiguration import Configuration

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
