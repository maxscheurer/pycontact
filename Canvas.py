'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''

from PyQt5.QtGui import (QColor, QPainter, QFont, QPixmap)
from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import QSize
import numpy as np

from biochemistry import *
from filters import *
from LabelView import LabelView


class ColorScheme:
    custom, bbsc = range(2)


class Canvas(QWidget):
    def __init__(self):
        super(QWidget, self).__init__()

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
        self.range = [0, 0]
        self.rangeFilterActive = False
        self.showHbondScores = False

    def paintEvent(self, event):

        qp = QPainter()
        qp.begin(self)

        # render pixmap to resolve performance issues
        if self.rendered:
            self.drawRenderedContact(event, qp)
        elif self.rendered is False and self.contacts:
            self.renderContact(False)
            self.rendered = True

        self.setMinimumSize(QSize(self.sizeX, self.sizeY))

        qp.end()

    def renderContact(self, generator):
        # startx = 90
        # orig_startx = startx
        start_text = 10
        rowheight = 22
        textoffset = 12
        blackColor = QColor(0, 0, 0)
        whiteColor = QColor(255, 255, 255)

        merge = self.merge
        offset = 10

        self.labelView.clean()
        self.labelView = LabelView(self.contacts)
        self.labelView.setParent(self)
        self.labelView.nsPerFrame = self.nsPerFrame
        self.labelView.threshold = self.threshold
        self.labelView.show()
        # startx has to be set according to maximum button length in labelview
        startx = np.max(self.labelView.buttonWidths) + 15
        orig_startx = startx

        # self.sizeX = (len(self.contacts[0].scoreArray) + startx) * offset
        # self.sizeY = len(self.contacts) * rowheight
        if self.rangeFilterActive:
            self.sizeX = startx + (len(self.contacts[0].scoreArray) + merge * 2) * offset / merge
        else:
            self.sizeX = startx + (len(
                self.contacts[0].scoreArray[self.range[0]:self.range[1]]) + merge * 2) * offset / merge

        self.sizeY = len(self.contacts) * rowheight

        self.pixmap = QPixmap(QSize(self.sizeX, self.sizeY))
        p = QPainter()

        if generator:
            p.begin(generator)
        else:
            p.begin(self.pixmap)

        p.fillRect(0, 0, self.sizeX, self.sizeY, whiteColor)

        row = 0
        self.alphaFactor = 50
        for c in self.contacts:
            bbScColor = BackboneSidechainContactType.colors[c.determineBackboneSidechainType()]
            i = 0
            if not self.showHbondScores:
                if self.rangeFilterActive:
                    rangedScores = c.scoreArray
                else:
                    rangedScores = c.scoreArray[self.range[0]:self.range[1]]
            else:
                hbarray = c.hbondFramesScan()
                self.alphaFactor = 100
                if self.rangeFilterActive:
                    rangedScores = hbarray
                else:
                    rangedScores = hbarray[self.range[0]:self.range[1]]
            while i < len(rangedScores):
                p.setPen(blackColor)
                merged_score = 0
                for j in range(merge):
                    if (i + j) >= len(rangedScores):
                        break
                    x = rangedScores[i + j]
                    merged_score += x
                merged_score = merged_score / merge
                alpha = merged_score * self.alphaFactor
                if alpha > 255:
                    alpha = 255
                if self.colorScheme == ColorScheme.bbsc:
                    p.setBrush(QColor(bbScColor[0], bbScColor[1], bbScColor[2], alpha))
                elif self.colorScheme == ColorScheme.custom:
                    color = QColor(self.customColor)
                    color.setAlpha(alpha)
                    p.setBrush(color)
                p.drawRect(startx, row, offset, 20)
                startx += offset
                i += merge
            startx = orig_startx
            row += rowheight

        if generator:
            row = 0
            for c in self.contacts:
                p.setPen(0)
                # print(ContactType.colors[c.contactType])
                p.setFont(QFont('Arial', 9))
                # string = c.resA + c.residA + "-" + c.resB + c.residB
                string = c.title
                p.setBrush(ContactType.qcolors[c.determine_ctype()])
                p.drawRect(0, row+3, orig_startx, rowheight-10)
                p.setPen(1)
                p.drawText(start_text, row + textoffset, string)
                row += rowheight

        p.end()

    def drawRenderedContact(self, event, qp):
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)
