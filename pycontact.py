import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton,
                             QFrame, QApplication)
import shelve
import numpy as np


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):

        centralWidget = QWidget()

        menubar = QMenuBar()
        menubar.addMenu('File')

        grid = QGridLayout()
        grid.setSpacing(10)

        self.painter = Canvas()
        grid.addWidget(self.painter, 1, 0)

        centralWidget.setLayout(grid)
        self.setCentralWidget(centralWidget)

        openFile = QAction(QIcon("open.png"), "Open", self)
        openFile.setShortcut("Ctrl+O")
        openFile.triggered.connect(self.showOpenFile)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu("&File")
        fileMenu.addAction(openFile)

        # contact = Contact("A","B",1,2,[1,2])
        # print(contact.type)
        # print(contact.framenumber())

        # User Defaults

        # self.defaults = shelve.open("defaultState")
        # if len(self.defaults) > 0:
        #     self.userField.setText(self.defaults["user"])
        #     self.user = self.defaults["user"]
        #     self.hostField.setText(self.defaults["host"])
        #     self.host = self.defaults["host"]
        #     self.passwdField.setText(self.defaults["passwd"])
        #     self.passwd = self.defaults["passwd"]
        #     self.cbSaveLogin.toggle()

        self.resize(1000, 550)
        self.center()
        self.setWindowTitle("pyContact")
        self.show()

    def showOpenFile(self):
        fnames = QFileDialog.getOpenFileNames(self, "Open file")
        for file in fnames[0]:
            self.file = file
            break
        lines = []
        print(self.file)
        with open(self.file, "r") as f:
            for line in f.readlines():
                lines.append(line)
        self.makeContactFromLines(lines)

    def makeContactFromLines(self, lines):
        self.contacts = []
        for l in lines:
            single_elements = l.split()
            resA = single_elements[0]
            residA = single_elements[1]
            resB = single_elements[2]
            residB = single_elements[3]
            frame_info = np.array(single_elements[4:], dtype=float)
            newContact = Contact(resA, residA, resB, residB, frame_info)
            self.contacts.append(newContact)
        print("new contacts: " + str(len(self.contacts)))
        self.painter.contacts = self.contacts
        self.painter.update()
    # self.painter.take_screenshot("test")



    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def pushLogin(self):

        self.host = self.hostField.text()
        self.user = self.userField.text()
        self.passwd = self.passwdField.text()

        try:
            self.ssh = paramiko.SSHClient()
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.ssh.connect(self.host, username=self.user, password=self.passwd)
            print("Successfully connected to " + self.host)

        except paramiko.AuthenticationException:
            print("Authentication error while connecting to " + host)

        if self.cbSaveLogin.isChecked():
            self.defaults = shelve.open("defaultState")
            self.defaults["user"] = self.user
            self.defaults["host"] = self.host
            self.defaults["passwd"] = self.passwd
            self.defaults.close()


class Canvas(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        # self.col = QColor(0, 0, 0)

        # self.square = QFrame(self)
        # self.square.setGeometry(0, 0, 1, 100)
        # self.square.setStyleSheet("QWidget { background-color: %s }" %
        # self.col.name())
        self.contacts = 0

    def paintEvent(self, event):

        qp = QPainter()
        qp.begin(self)
        self.drawContact(event, qp)
        qp.end()

    def take_screenshot(self, filename):
        print(filename)
        p = QPixmap(self.size())
        self.render(p)
        p.save(filename, 'png', 100)

    def drawContact(self, event, qp):

        startx = 100
        orig_startx = startx
        start_text = 10
        textoffset = 13
        rowheight = 22
        offset = 10
        color = QColor(0, 0, 0)
        if self.contacts:
            row = 5
            for c in self.contacts:
                for x in c.scoreArray:
                    qp.setPen(color)
                    qp.setBrush(QColor(0, 200, 0, x * 50))
                    # qp.setBrush(QColor(((1-x)*50), x*50,0))
                    qp.drawRect(startx, row, offset, 20)
                    startx += (offset)
                    qp.setFont(QFont('Arial', 9))
                    string = c.resA + c.residA + "-" + c.resB + c.residB
                # qp.drawText(start_text, row+textoffset ,string)
                startx = orig_startx
                row += rowheight


class AnalysisView(QWidget):
    """docstring for AnalysisView"""

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.contacts = 0


class Contact:
    def __init__(self, resA, residA, resB, residB, scoreArray):
        self.resA = resA
        self.resB = resB
        self.residA = residA
        self.residB = residB
        self.scoreArray = scoreArray
        self.type = determine_ctype(self.resA, self.resB)

    def framenumber(self):
        return len(self.scoreArray)

def determine_ctype(resA, resB):
    return 0


def main():
    app = QApplication(sys.argv)
    mainWindow = MainWindow()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
