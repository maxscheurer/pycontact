import sip

from PyQt5.QtGui import (QColor, QPainter, QFont, QPixmap, QPaintEvent)
from PyQt5.QtWidgets import QWidget, QLabel, QSizePolicy
from PyQt5.QtCore import QSize, QRect
from PyQt5.QtCore import pyqtSignal, QObject
import numpy as np

from ..core.ContactFilters import *
from .LabelView import LabelView


class StretchedLabel(QLabel):
    def __init__(self, *args, **kwargs):
        QLabel.__init__(self, *args, **kwargs)
        self.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)

    def scaleLabel(self, factor):
        font = self.font()
        font.setPixelSize(self.height() * factor)
        self.setFont(font)


class ColorScheme:
    """Enum the color scheme, either custom or backbone-sidechain type."""
    custom, bbsc = range(2)


class TrackCanvas(QWidget, QObject):
    """Canvas where track molecule analysis results are drawn."""

    def __init__(self):
        super(QWidget, self).__init__()
        self.rendered = 0
        self.contacts = None
        self.labels = []
        self.maximalContactsPerRow = 1
        self.sizeX = 0
        self.sizeY = 0

    def draw_labels(self):
        if self.contacts is not None:
            self.clean()
            self.renderContact(False)

    def paintEvent(self, event):
        qp = QPainter()
        qp.begin(self)
        # render pixmap to resolve performance issues
        # if self.rendered:
        #     self.drawRenderedContact(qp)
        # elif self.rendered is False and self.contacts:
        #     print("render contact if")
        #     self.renderContact(False)
        # self.rendered = True

        self.setMinimumSize(QSize(self.sizeX, self.sizeY))

        qp.end()


    def clean(self):
        """Delete all labels."""
        allLabels = self.findChildren(StretchedLabel)
        for child in allLabels:
            sip.delete(child)


    def renderContact(self, generator):
        """Render the tracking timeline."""

        blackColor = QColor(0, 0, 0)
        whiteColor = QColor(255, 255, 255)

        # self.pixmap = QPixmap(QSize(self.sizeX, self.sizeY))
        #p = QPainter()

        # if generator:
        #     p.begin(generator)
        # else:
        #     p.begin(self.pixmap)

        # p.fillRect(0, 0, self.sizeX, self.sizeY, whiteColor)
        #p.end()
        row = 0
        rownumber = 0
        colnumber = 0
        offsetY = 10
        offsetX = 50
        colwidth = 80
        colnumberMax = 10
        lheight = 20
        totalWidth = offsetX
        totalHeight = offsetY
        rowheight = 0

        frameLabels = []
        # print("attempt to paint", len(self.contacts))
        frameCounter = 1
        for frame in self.contacts:
            # print("painting labels", len(frame))
            currentLabels = []
            lnumber = 0
            # print(frame)
            if len(frame) == 0:
                currentLabels.append(StretchedLabel("empty"))
                currentLabels[-1].setParent(self)
                width = currentLabels[-1].width()
                # print("width", width)
                currentLabels[-1].move(offsetX + colnumber*colwidth, offsetY + lnumber*lheight + totalHeight)
                currentLabels[-1].setFont(QFont('Arial', 10))
                currentLabels[-1].show()
                rowheight = totalHeight + self.maximalContactsPerRow*lheight

            for label in frame:
                stylesheet = "border: 0px solid #222222; background-color: rgba(255,0,0,%d);" % int(label[1]*2)
                currentLabels.append(StretchedLabel(label[0]))
                currentLabels[-1].setAutoFillBackground(True)
                currentLabels[-1].setStyleSheet(stylesheet)
                currentLabels[-1].setParent(self)
                width = currentLabels[-1].width()
                # print("width", width)
                currentLabels[-1].move(offsetX + colnumber*colwidth, offsetY + lnumber*lheight + totalHeight)
                currentLabels[-1].setFont(QFont('Arial', 10))
                currentLabels[-1].show()
                lnumber += 1
                if lnumber >= self.maximalContactsPerRow:
                    break
            rowheight = totalHeight + self.maximalContactsPerRow*lheight


            self.labels.append(currentLabels)
            colnumber += 1

            if colnumber >= colnumberMax:
                frameCounter += colnumber
                colnumber = 0
                rownumber += 1
                offsetY += 15
                totalHeight = rowheight
                if frameCounter < len(self.contacts):
                    frameLabels.append(StretchedLabel("%d" % frameCounter))
                    frameLabels[-1].setParent(self)
                    width = frameLabels[-1].width()
                    # print("width", width)
                    frameLabels[-1].move(10 + colnumber*colwidth, offsetY + 0.4*self.maximalContactsPerRow*lheight + totalHeight)
                    frameLabels[-1].setFont(QFont('Arial', 10))
                    frameLabels[-1].show()

        self.sizeX = colnumberMax * colwidth + offsetX
        self.sizeY = offsetY + rownumber * self.maximalContactsPerRow * lheight

    def drawRenderedContact(self, qp):
        """Draws the rendered contact to the canvas."""
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)
