import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QVBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton,
                             QFrame, QApplication, QSizePolicy)
import shelve
import numpy as np
import gui
from functools import partial


class MainWindow(QMainWindow, gui.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)

        self.setWindowTitle("pyContact")

        self.painter = Canvas()
        self.scrollArea.setWidget(self.painter)
        self.actionOpen.triggered.connect(self.pushOpen)

    def pushOpen(self):
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
        self.sizeX = 0
        self.sizeY = 0
        self.rendered = False
        self.pixmap = 0

    def paintEvent(self, event):

        qp = QPainter()
        qp.begin(self)

        #render pixmap to resolve performance issues
        if self.rendered:
            self.drawRenderedContact(event, qp)
        elif self.rendered == False and self.contacts:
            self.renderContact()
            self.rendered = True

        self.setMinimumSize(QSize(self.sizeX, self.sizeY))

        qp.end()

    def take_screenshot(self, filename):
        print(filename)
        p = QPixmap(self.size())
        self.render(p)
        p.save(filename, 'png', 100)

    def renderContact(self):
        startx = 80
        orig_startx = startx
        start_text = 10
        textoffset = 5
        rowheight = 22
        offset = 10
        blackColor = QColor(0, 0, 0)
        whiteColor = QColor(255, 255, 255)


        self.sizeX = (len(self.contacts[0].scoreArray) + startx) * offset
        self.sizeY = len(self.contacts) * rowheight

        self.pixmap = QPixmap(QSize(self.sizeX, self.sizeY))
        p = QPainter()
        p.begin(self.pixmap)
        p.fillRect(0, 0, self.sizeX, self.sizeY, whiteColor)

        row = 0
        for c in self.contacts:
            # print total time of contact in nanoseconds
            #print(c.total_time(1, 0.5))
            for x in c.scoreArray:
                p.setPen(blackColor)
                p.setBrush(QColor(0, 200, 0, x * 50))
                # qp.setBrush(QColor(((1-x)*50), x*50,0))
                p.drawRect(startx, row, offset, 20)
                startx += (offset)
                # p.setFont(QFont('Arial', 9))
                # p.drawText(start_text, row+textoffset ,string)
            startx = orig_startx
            row += rowheight

        p.end()
        self.pixmap.save("test", 'png', 100)
        self.analysis = AnalysisView(self.contacts)
        self.analysis.setParent(self)
        self.analysis.show()

    def drawRenderedContact(self, event, qp):
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)


class AnalysisView(QWidget):
    """docstring for AnalysisView"""

    def __init__(self, contacts):
        super().__init__()
        self.contacts = contacts
        self.initUI()

    def initUI(self):
        startx = 80
        orig_startx = startx
        start_text = 10
        textoffset = 5
        rowheight = 22
        row = 0
        self.buttons = []
        for c in self.contacts:
            cindex = self.contacts.index(c)
            string = c.resA +c.residA + "-" + c.resB + c.residB
            self.buttons.append(QPushButton(string))
            self.buttons[-1].setStyleSheet("border: 0px solid #222222")
            self.buttons[-1].clicked.connect(partial(self.handleButton, data=cindex))
            self.buttons[-1].setParent(self)
            self.buttons[-1].move(start_text, row + textoffset)
            self.buttons[-1].setFont(QFont('Arial', 9))
            self.buttons[-1].show()
            startx = orig_startx
            row += rowheight

    def handleButton(self, data):
        print('index clicked: '+ str(data))




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

    def total_time(self, ns_per_frame, threshold):
        time = 0
        for score in self.scoreArray:
            if score > threshold:
                time += ns_per_frame
        self.ttime = time
        return self.ttime

def determine_ctype(resA, resB):
    return 0


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()

class ContactType:
    hbond, saltbr, hydrophobic = range(3)

class ResidueType:
    positive, negative, unpolar, polar = range(4)



if __name__ == '__main__':
    main()
