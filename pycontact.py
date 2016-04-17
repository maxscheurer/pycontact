import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget, QDialog,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QVBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton,
                             QFrame, QApplication, QSizePolicy)
import numpy as np
import gui
from biochemistry import *
from inputreader import *
from functools import partial

# threshold: contact counted or not
threshold = 1
# nanoseconds per frame
nsperframe = 1

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
        self.contacts = makeContactFromLines(lines)
        print("new contacts: " + str(len(self.contacts)))
        self.painter.contacts = self.contacts
        self.painter.update()


class Canvas(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
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
            for x in c.scoreArray:
                p.setPen(blackColor)
                p.setBrush(QColor(0, 200, 0, x * 80))
                p.drawRect(startx, row, offset, 20)
                startx += (offset)
            startx = orig_startx
            row += rowheight

        p.end()
        self.pixmap.save("test", 'png', 100)
        self.labelView = LabelView(self.contacts)
        self.labelView.setParent(self)
        self.labelView.show()

    def drawRenderedContact(self, event, qp):
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)


class LabelView(QWidget):
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
            self.buttons.append(QPushButton(c.title))
            stylesheet = "border: 0px solid #222222; background-color: "+ ContactType.colors[c.type] + " ;"
            self.buttons[-1].setStyleSheet(stylesheet)
            self.buttons[-1].clicked.connect(partial(self.handleButton, data=cindex))
            self.buttons[-1].setParent(self)
            self.buttons[-1].move(start_text, row + textoffset)
            self.buttons[-1].setFont(QFont('Arial', 9))
            self.buttons[-1].show()
            startx = orig_startx
            row += rowheight

    def handleButton(self, data):
        print('index clicked: '+ str(data))
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)
        # b1 = QPushButton("ok", d)
        # b1.move(50, 50)
        contact = self.contacts[data]
        timeLabel = QLabel(str(contact.total_time(nsperframe,threshold)))
        timeTitleLabel = QLabel("total time [ns]:")
        grid.addWidget(timeTitleLabel, 0,0)
        grid.addWidget(timeLabel,0,1)
        d.setWindowTitle(contact.title)
        d.resize(200,100)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
