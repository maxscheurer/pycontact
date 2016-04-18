import sys, sip
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget, QDialog, QTabWidget,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QVBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton,
                             QFrame, QApplication, QSizePolicy)
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import gui
from settings import *
from biochemistry import *
from inputreader import *
from filters import *
from functools import partial

class MainWindow(QMainWindow, gui.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.contacts = []
        self.setupUi(self)

        self.setWindowTitle("pyContact")
        self.mergeSlider.valueChanged.connect(self.mergeValueChanged)

        self.painter = Canvas()
        self.scrollArea.setWidget(self.painter)
        self.actionOpen.triggered.connect(self.pushOpen)

        self.settingsView = SettingsTabWidget()
        self.settingsView.applySettingsButton.clicked.connect(self.updateSettings)
        self.settingsView.applyFilterButton.clicked.connect(self.updateFilters)

        self.alphaSlider.setValue(50)
        self.alphaSlider.valueChanged.connect(self.alphaValueChanged)

        self.statisticsButton.clicked.connect(self.showStatistics)

        self.updateSettings()
        self.updateFilters()
        self.openPreferencesButton.clicked.connect(self.openPrefs)

    def updateSettings(self):
        self.painter.nsPerFrame = float(self.settingsView.nsPerFrameField.text())
        self.painter.threshold = float(self.settingsView.thresholdField.text())
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def updateFilters(self):
        print("filter update")
        # total time filter
        totalTimeActive = self.settingsView.activeTotalTimeCheckbox.isChecked()
        scoreActive = self.settingsView.activeScoreCheckbox.isChecked()
        filterActive = (totalTimeActive or scoreActive)
        if  filterActive:
            print("filter act.")
            if len(self.contacts) > 0:
                filteredContacts = self.contacts
                if totalTimeActive:
                    operator = self.settingsView.compareTotalTimeDropdown.currentText()
                    value = float(self.settingsView.totalTimeField.text())
                    filter = TotalTimeFilter("tottime", operator, value)
                    filteredContacts = filter.filterContacts(filteredContacts)
                if scoreActive:
                    operator = self.settingsView.compareScoreDropdown.currentText()
                    value = float(self.settingsView.scoreField.text())
                    filter = ScoreFilter("score", operator, value)
                    filteredContacts = filter.filterContacts(filteredContacts)
                self.painter.contacts = filteredContacts
                self.painter.rendered = False
                self.painter.update()
                self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))
                if len(filteredContacts) == 0:
                    self.painter.labelView.clean()
        else:
            self.painter.contacts = self.contacts
            self.painter.rendered = False
            self.painter.update()
            self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))


    def openPrefs(self):
        self.settingsView.show()

    def showStatistics(self):
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)

        numberTitleLabel = QLabel("total number of contacts:")
        numberLabel = QLabel(str(len(self.contacts)))

        numberFramesTitleLabel = QLabel("number of frames:")
        numberFramesLabel = QLabel(str(len(self.contacts[0].scoreArray)))

        meanTitleLabel = QLabel("mean contact score:")
        meanLabel = QLabel(str(mean_score_of_contactArray(self.contacts)))

        grid.addWidget(numberTitleLabel, 0, 0)
        grid.addWidget(numberLabel, 0, 1)
        grid.addWidget(numberFramesTitleLabel, 1, 0)
        grid.addWidget(numberFramesLabel, 1, 1)
        grid.addWidget(meanTitleLabel, 2, 0)
        grid.addWidget(meanLabel, 2, 1)
        # contactPlot = ContactPlotter(None, width=4, height=2, dpi=80)
        # contactPlot.plot_contact_figure(contact)
        # grid.addWidget(contactPlot, 2, 0, 1, 2)
        d.setWindowTitle("Statistics")
        d.resize(600, 450)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

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
        # self.painter.update()
        # set max slider value to frame number!
        self.mergeSlider.setMaximum(len(self.contacts[0].scoreArray)/12)
        self.updateSettings()
        self.updateFilters()

    def mergeValueChanged(self):
        self.painter.merge = self.mergeSlider.value()
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

    def alphaValueChanged(self):
        self.painter.alphaFactor = self.alphaSlider.value()
        self.painter.rendered = False
        self.painter.update()
        self.painter.paintEvent(QPaintEvent(QRect(0, 0, self.painter.sizeX, self.painter.sizeY)))

class SettingsTabWidget(QTabWidget, Ui_settingsWindowWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QTabWidget, self).__init__(parent)
        self.setupUi(self)


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
        self.merge = 1
        self.labelView = LabelView([])
        self.alphaFactor = 50
        self.contacts = []
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
        blackColor = QColor(0, 0, 0)
        whiteColor = QColor(255, 255, 255)

        merge = self.merge
        offset = 10 / merge

        self.sizeX = (len(self.contacts[0].scoreArray) + startx) * offset
        self.sizeY = len(self.contacts) * rowheight

        self.pixmap = QPixmap(QSize(self.sizeX, self.sizeY))
        p = QPainter()
        p.begin(self.pixmap)
        p.fillRect(0, 0, self.sizeX, self.sizeY, whiteColor)

        row = 0
        for c in self.contacts:
            i = 0
            while i < len(c.scoreArray):
                p.setPen(blackColor)
                merged_score = 0
                for j in range(merge):
                    if (i+j) >= len(c.scoreArray):
                        break
                    x = c.scoreArray[i+j]
                    merged_score += x
                merged_score = merged_score / merge
                alpha = merged_score * self.alphaFactor
                if alpha > 255:
                    alpha = 255
                p.setBrush(QColor(0, 200, 0, alpha))
                p.drawRect(startx, row, offset*merge, 20)
                startx += (offset*merge)
                i += merge
            startx = orig_startx
            row += rowheight

        p.end()
        # self.pixmap.save("test", 'png', 100)
        self.labelView.clean()
        self.labelView = LabelView(self.contacts)
        self.labelView.setParent(self)
        self.labelView.nsPerFrame = self.nsPerFrame
        self.labelView.threshold = self.threshold
        self.labelView.show()

    def drawRenderedContact(self, event, qp):
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)


class LabelView(QWidget):
    """docstring for AnalysisView"""

    def __init__(self, contacts):
        super().__init__()
        self.contacts = contacts
        self.initUI()

    def clean(self):
        allLabels = self.findChildren(QPushButton)
        for child in allLabels:
            sip.delete(child)

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
        timeLabel = QLabel(str(contact.total_time(self.nsPerFrame,self.threshold)))
        thresholdLabel = QLabel(str(self.threshold))
        timeTitleLabel = QLabel("total time [ns]:")
        thresholdTitleLabel = QLabel("current threshold:")

        meanTitleLabel = QLabel("mean score:")
        meanLabel = QLabel(str(contact.mean_score()))
        grid.addWidget(timeTitleLabel, 0,0)
        grid.addWidget(timeLabel,0,1)
        grid.addWidget(thresholdTitleLabel,1,0)
        grid.addWidget(thresholdLabel, 1, 1)
        grid.addWidget(meanTitleLabel, 2, 0)
        grid.addWidget(meanLabel, 2, 1)
        contactPlot = ContactPlotter(None, width=4, height=2, dpi=80)
        contactPlot.plot_contact_figure(contact)
        grid.addWidget(contactPlot,3,0,1,2)
        d.setWindowTitle(contact.title)
        d.resize(600,450)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()


class MplPlotter(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class ContactPlotter(MplPlotter):
    """Simple canvas with a sine plot."""

    def plot_contact_figure(self, contact):
        self.axes.plot(contact.scoreArray)
        self.axes.set_xlabel("frame")
        self.axes.set_ylabel("score")


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()
