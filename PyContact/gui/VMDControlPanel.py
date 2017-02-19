from PyQt5.QtWidgets import QWidget, QGridLayout, QLabel, QPushButton, QComboBox, QLineEdit, QCheckBox
from PyQt5 import QtCore

from ..cy_modules import wrap_vmd as vmd

class VMDControlPanel(QWidget):
    """docstring for VMDControlPanel"""

    def __init__(self):
        super(QWidget,self).__init__()
        self.initUI()


    def initUI(self):
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.setWindowTitle("VMD Control Panel")
        self.resize(640, 444)
        self.grid.setGeometry(QtCore.QRect(10, 10, 621, 431))

        self.startButton = QPushButton("Start VMD")
        self.startButton.clicked.connect(self.pushStartVMD)
        self.grid.addWidget(self.startButton, 0, 0)
        self.startButton.setEnabled(True)

        self.stopButton = QPushButton("Stop VMD")
        self.stopButton.clicked.connect(self.pushStopVMD)
        self.grid.addWidget(self.stopButton, 0, 1)
        self.stopButton.setEnabled(False)

# just for testing purposes
        self.commandButton = QPushButton("Send command")
        self.commandButton.clicked.connect(self.sendCommand)
        self.grid.addWidget(self.commandButton, 2, 0)
        self.commandButton.setEnabled(False)

        self.commandField = QLineEdit()
        self.grid.addWidget(self.commandField, 1, 0, 1, 2)

    def pushStartVMD(self):
        self.startButton.setEnabled(False)
        self.stopButton.setEnabled(True)
        self.commandButton.setEnabled(True)
        vmd.start(5050,"vmd")

    def pushStopVMD(self):
        self.startButton.setEnabled(True)
        self.stopButton.setEnabled(False)
        self.commandButton.setEnabled(False)
        vmd.stop()

    def sendCommand(self):
        vmd.send_command(self.commandField.text().rstrip())
