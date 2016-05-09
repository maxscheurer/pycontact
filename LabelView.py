import sys, sip, copy
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget, QDialog, QTabWidget, QButtonGroup,
                             QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
                             QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QVBoxLayout)

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtGui import (QColor, QPainter, QFont)
from PyQt5.QtWidgets import (QWidget, QPushButton, QRadioButton,
                             QFrame, QApplication, QSizePolicy)
from PyQt5.QtSvg import QSvgGenerator
import numpy as np
from settings import *
from biochemistry import *
from functools import partial
from Plotters import *

class LabelView(QWidget):
    """docstring for AnalysisView"""

    def __init__(self, contacts):
        super(QWidget,self).__init__()
        self.contacts = contacts
        self.initUI()

    def clean(self):
        allLabels = self.findChildren(QPushButton)
        for child in allLabels:
            sip.delete(child)

    def initUI(self):
        start_text = 10
        textoffset = 5
        rowheight = 22
        row = 0
        self.buttons = []
        for c in self.contacts:
            cindex = self.contacts.index(c)
            # TODO: implement title in AccumulatedContact
            self.buttons.append(QPushButton(c.title))
            # TODO: implement ContactType in AccumulatedContact
            # stylesheet = "border: 0px solid #222222; background-color: " + ContactType.colors[c.contactType] + " ;"
            stylesheet = "border: 0px solid #222222; background-color: " + ContactType.colors[0] + " ;"
            self.buttons[-1].setStyleSheet(stylesheet)
            self.buttons[-1].clicked.connect(partial(self.handleButton, data=cindex))
            self.buttons[-1].setParent(self)
            self.buttons[-1].move(start_text, row + textoffset)
            self.buttons[-1].setFont(QFont('Arial', 9))
            self.buttons[-1].show()
            row += rowheight

    def handleButton(self, data):
        # print('index clicked: '+ str(data))
        d = QDialog()
        grid = QGridLayout()
        d.setLayout(grid)
        contact = self.contacts[data]
        timeLabel = QLabel(str(contact.total_time(self.nsPerFrame, self.threshold)))
        thresholdLabel = QLabel(str(self.threshold))
        timeTitleLabel = QLabel("total time [ns]:")
        thresholdTitleLabel = QLabel("current threshold:")

        backboneSidechainTitleLabel = QLabel("bb/sc score (A)")
        backboneSidechainTitleLabel2 = QLabel("bb/sc score (B)")

        backboneSidechainLabel = QLabel("%.2f/%.2f" % (contact.bb1, contact.sc1))
        backboneSidechainLabel2 = QLabel("%.2f/%.2f" % (contact.bb2, contact.sc2))

        meanLifeTimeTitleLabel = QLabel("mean lifetime:")
        meanLifeTimeLabel = QLabel(str(contact.mean_life_time(self.nsPerFrame, self.threshold)))

        medianLifeTimeTitleLabel = QLabel("median lifetime:")
        medianLifeTimeLabel = QLabel(str(contact.median_life_time(self.nsPerFrame, self.threshold)))

        meanTitleLabel = QLabel("mean score:")
        meanLabel = QLabel(str(contact.mean_score()))

        medianTitleLabel = QLabel("median score:")
        medianLabel = QLabel(str(contact.median_score()))

        grid.addWidget(timeTitleLabel, 0, 0)
        grid.addWidget(timeLabel, 0, 1)
        grid.addWidget(thresholdTitleLabel, 1, 0)
        grid.addWidget(thresholdLabel, 1, 1)
        grid.addWidget(meanTitleLabel, 2, 0)
        grid.addWidget(meanLabel, 2, 1)

        grid.addWidget(medianTitleLabel, 2, 2)
        grid.addWidget(medianLabel, 2, 3)

        grid.addWidget(backboneSidechainTitleLabel, 3, 0)
        grid.addWidget(backboneSidechainLabel, 3, 1)
        grid.addWidget(backboneSidechainTitleLabel2, 4, 0)
        grid.addWidget(backboneSidechainLabel2, 4, 1)

        grid.addWidget(meanLifeTimeTitleLabel, 5, 0)
        grid.addWidget(meanLifeTimeLabel, 5, 1)
        grid.addWidget(medianLifeTimeTitleLabel, 5, 2)
        grid.addWidget(medianLifeTimeLabel, 5, 3)

        contactPlot = ContactPlotter(None, width=4, height=2, dpi=80)
        contactPlot.plot_contact_figure(contact)
        grid.addWidget(contactPlot, 6, 0, 1, 4)
        d.setWindowTitle(contact.title)
        d.resize(650, 700)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()