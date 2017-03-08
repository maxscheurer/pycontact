from socket import *
import subprocess
from pkg_resources import resource_filename as res

from PyQt5.QtWidgets import QWidget, QGridLayout, QLabel, QPushButton, QLineEdit
from PyQt5 import QtCore

from Dialogues import TopoTrajLoaderDialog
from ..core.Biochemistry import *

# from ..cy_modules import wrap_vmd as vmd


class VMDCommands():

    @staticmethod
    def translateSelections(mdanalysis_text):
        txt = mdanalysis_text.replace("segid","segname")
        txt = txt.replace("-"," to ")
        #around x y--> within x of y
        lst = txt.split(" ")
        idcs = [i for i, j in enumerate(lst) if j == 'around']
        for idx in idcs:
            lst.insert(idx+2,"of")
        txt = " ".join(lst)
        txt = txt.replace("around","within")
        return txt


    @staticmethod
    def gotoFrame(frame):
        return "animate goto %s" % str(frame)

    @staticmethod
    def styleBackbone():
        return """
        mol addrep 0
        mol modstyle top 0 NewCartoon 0.300000 10.000000 4.100000 0
        mol modselect top 0 backbone
        mol modcolor top 0 Chain
        """

    @staticmethod
    def addSelection(sel, representations, colorID):
        idx = str(len(representations))
        representations.append(sel)
        return ["""
        mol addrep top
        mol modstyle %s top Licorice
        mol modselect %s top (%s)
        mol modcolor %s top ColorID %d
        """ % (idx, idx, sel, idx, colorID), representations]

    @staticmethod
    def removeReps(index):
        # delrep rep_number molecule_number
        return """
        mol delrep %s top
        """ % index

    @staticmethod
    def resetView():
        return "display resetview"


class VMDTcp():
    """Interface to VMD via tcp"""
    def __init__(self):
        self.commands = VMDCommands()
        self.rctl = res(__name__, './remote_ctl.tcl')
        self.HOST = 'localhost'
        self.PORT = 5050
        self.ADDR = (self.HOST, self.PORT)

    def attemptConnection(self):
        self.tcpClientSocket = socket(AF_INET, SOCK_STREAM)
        self.tcpClientSocket.connect(self.ADDR)

    def start(self):
        subprocess.Popen(["vmd","-e",self.rctl])
        try:
            self.attemptConnection()
        except:
            return -1

    def send_command(self,cmd):
        self.tcpClientSocket.send(str(cmd+"\n"))

    def stop(self):
        self.send_command("quit")
        self.tcpClientSocket.close()

class VMDControlPanel(QWidget):
    """docstring for VMDControlPanel"""

    def __init__(self):
        super(QWidget,self).__init__()
        self.initUI()
        self.representations = []
        self.connected = False


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

        self.loadTopoTrajButton = QPushButton("Load Molecules into VMD")
        self.loadTopoTrajButton.clicked.connect(self.loadTopoTraj)
        self.grid.addWidget(self.loadTopoTrajButton, 2, 1)
        self.loadTopoTrajButton.setEnabled(False)

        self.connectButton = QPushButton("Connect to VMD")
        self.connectButton.clicked.connect(self.pushConnectVMD)
        self.grid.addWidget(self.connectButton, 0, 1)
        self.connectButton.setEnabled(False)

        self.infoLabel = QLabel("No VMD session running")
        self.grid.addWidget(self.infoLabel, 3, 0, 2, 2)

        self.stopButton = QPushButton("Stop VMD")
        self.stopButton.clicked.connect(self.pushStopVMD)
        self.grid.addWidget(self.stopButton, 0, 2)
        self.stopButton.setEnabled(False)

# just for testing purposes
        self.commandButton = QPushButton("Send command")
        # self.commandButton.clicked.connect(self.sendCommand)
        # self.grid.addWidget(self.commandButton, 2, 0)
        # self.commandButton.setEnabled(False)

        self.commandField = QLineEdit()
        self.grid.addWidget(self.commandField, 1, 0, 1, 2)
        self.vmd = VMDTcp()

    def addRep(self, txt):
        self.representations.append(txt)

    def prepareVMDWithTopoTraj(self,top,traj):
        self.vmd.send_command("mol new %s"  % top)
        self.vmd.send_command(self.vmd.commands.removeReps(0))
        self.vmd.send_command("mol addfile %s waitfor all" % traj)
        self.vmd.send_command(self.vmd.commands.styleBackbone())
        self.vmd.send_command(self.vmd.commands.resetView())
        self.vmd.send_command(self.vmd.commands.gotoFrame(0))
        self.addRep("initialBB")

    def loadTopoTraj(self):
        topoloader = TopoTrajLoaderDialog()
        cfg, result = topoloader.getConfig()
        if result == 1:
            self.prepareVMDWithTopoTraj(cfg[0], cfg[1])

    def gotoVMDFrame(self, frame):
        self.vmd.send_command(self.vmd.commands.gotoFrame(frame))

    def updateSelections(self, main_sel1, main_sel2, cont_list):
        for s in reversed(range(1, len(self.representations) + 1)):
            self.vmd.send_command(self.vmd.commands.removeReps(s))
        self.representations = self.representations[:1]

        sel1 = self.vmd.commands.translateSelections(main_sel1) + " and ("
        sel2 = self.vmd.commands.translateSelections(main_sel2) + " and ("
        for cont in cont_list:
            currentSel1 = []
            index = 0
            for item in cont.key1:
                if item != "none":
                    currentSel1.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                index += 1
            currentSel1String = " and ".join(currentSel1)
            sel1 += currentSel1String + " or "
            currentSel2 = []
            index = 0
            for item in cont.key2:
                if item != "none":
                    currentSel2.append(AccumulationMapIndex.vmdsel[index] + " " + item)
                index += 1
            currentSel2String = " and ".join(currentSel2)
            sel2 += currentSel2String + " or "
        sel1 = sel1[:-3] + ")"
        sel2 = sel2[:-3] + ")"

        # print(self.vmd.commands.addSelection(sel1, 3))
        sel1command, self.representations = self.vmd.commands.addSelection(sel1, self.representations, 3)
        self.vmd.send_command(sel1command)
        sel2command, self.representations = self.vmd.commands.addSelection(sel2, self.representations, 4)
        self.vmd.send_command(sel2command)
        # self.vmd.send_command(self.vmd.commands.styleBackbone())
        # self.vmd.send_command(self.vmd.commands.resetView())

    def updateInfoLabel(self, txt):
        self.infoLabel.setText(txt)

    def pushConnectVMD(self):
        try:
            self.vmd.attemptConnection()
            self.updateInfoLabel("Connection established")
            self.loadTopoTrajButton.setEnabled(True)
            self.connected = True
        except Exception as e:
             self.updateInfoLabel("Could not connect to VMD!\nTry to connect using the Connect button\n when VMD is opened.")

    def pushStartVMD(self):
        self.startButton.setEnabled(False)
        self.stopButton.setEnabled(True)
        self.commandButton.setEnabled(True)
        res = self.vmd.start()
        if res == -1:
            self.connectButton.setEnabled(True)
            self.updateInfoLabel("Could not connect to VMD!\nTry to connect using the Connect button\n when VMD is opened.")
        else:
            self.loadTopoTrajButton.setEnabled(True)
            self.connected = True


    def pushStopVMD(self):
        self.representations = []
        self.updateInfoLabel("VMD stopped.")
        self.connectButton.setEnabled(False)
        self.startButton.setEnabled(True)
        self.stopButton.setEnabled(False)
        self.commandButton.setEnabled(False)
        self.loadTopoTrajButton.setEnabled(False)
        self.vmd.stop()
        self.connected = False

    def sendCommand(self):
        self.vmd.send_command(self.commandField.text())
        self.commandField.setText("")
