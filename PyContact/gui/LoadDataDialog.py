from PyQt5.QtWidgets import QDialog, QLineEdit, QFileDialog
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QVariant
from PyQt5.QtCore import Qt

from .LoadData import Ui_Dialog
from ..core.LoadConfiguration import Configuration
from ..exampleData.datafiles import DCD, PSF

class LoadDataDialog(QDialog, Ui_Dialog):
    """docstring for [object Object]."""
    def __init__(self, parent=None):
        super(LoadDataDialog, self).__init__(parent)
        self.setupUi(self)
        production = 0


        self.cutoffField = QLineEdit("5.0")
        posDoubleValidator = QDoubleValidator()
        posDoubleValidator.setBottom(0)
        self.cutoffField.setValidator(posDoubleValidator)
        self.cutoffAngleField = QLineEdit("120")
        self.cutoffAngleField.setValidator(posDoubleValidator)
        self.cutoffHbondField = QLineEdit("2.5")
        self.cutoffHbondField.setValidator(posDoubleValidator)

        self.trajectories = []
        self.topology = ""

        self.tableModel = TrajectoryTableModel()
        self.trajectoryTableView.setModel(self.tableModel)

        self.buttonTopology.clicked.connect(self.pick_topology)
        self.buttonTrajectory.clicked.connect(self.pick_trajectory)

        if not production:
            self.selection1Field.setText("segid RN11")
            self.selection2Field.setText("segid UBQ")

    def pick_topology(self):
        """Pick topology file."""
        topology = QFileDialog.getOpenFileNames(self, "Open topology")
        for file in topology[0]:
            self.topology = file
            self.topFileDisplay.setText(self.topology.split("/")[-1])
            break

    def pick_trajectory(self):
        """Pick trajectory file."""
        trajs = QFileDialog.getOpenFileNames(self, "Open trajectory")
        for file in trajs[0]:
            self.trajectories.append(file)
        self.tableModel.update(self.trajectories)
        self.trajectoryTableView.resizeColumnsToContents()

    def configuration(self):
        """Returns the chosen configuration."""
        config = Configuration(self.topology, self.trajectories,
                               float(self.cutoffField.text()),
                               float(self.cutoffHbondField.text()),
                               float(self.cutoffAngleField.text()),
                               self.selection1Field.text(),
                               self.selection2Field.text())
        return config

    @staticmethod
    def getConfig(parent=None):
        """Static method to create the dialog and return (date, time, accepted)."""
        dialog = LoadDataDialog(parent)
        result = dialog.exec_()
        config = dialog.configuration()
        return config, result == QDialog.Accepted


class TrajectoryTableModel(QAbstractTableModel):
    def __init__(self, parent=None, *args):
        super(TrajectoryTableModel, self).__init__()
        self.trajs = []

    def update(self, trajs):
        self.layoutAboutToBeChanged.emit()
        self.trajs = trajs
        self.layoutChanged.emit()

    def rowCount(self, parent):
        return len(self.trajs)

    def columnCount(self, parent):
        return 1

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self.trajs[index.row()].split("/")[-1]
        else:
            return None

    def flags(self, index):
        return Qt.ItemIsEnabled
