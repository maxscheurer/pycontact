from PyQt5.QtGui import (QColor, QPainter, QFont, QPixmap, QPaintEvent)
from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import QSize, QRect
from PyQt5.QtCore import pyqtSignal, QObject
import numpy as np

from ..core.ContactFilters import *
from .LabelView import LabelView


class ColorScheme:
    """Enum the color scheme, either custom or backbone-sidechain type."""
    custom, bbsc = range(2)


class Canvas(QWidget, QObject):
    """Canvas where contact analysis results are drawn."""
    clickedRowSignal = pyqtSignal()
    clickedColumnSignal = pyqtSignal()

    def __init__(self):
        super(QWidget, self).__init__()

        self.clickedRow = 0
        self.clickedColumn = 0
        self.rendered = False
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
        self.vismode = False
        self.timeLineXOrigin = 0
        self.clickedRow = -1
        self.clickedColumn = -1
        self.offset = -1
        self.globalClickedRow = -1
        self.timeLineXOrigin = 0
        self.rowh = 1
        self.endOfTimeLine = 0

    def mousePressEvent(self, event):
        pos = event.pos()
        x, y = pos.x(), pos.y()
        self.clickedRow = -1
        if self.vismode and self.timeLineXOrigin < x < self.endOfTimeLine:
            self.clickedRow = int(y / self.rowh - 1)  # -1 because of frame number line
            self.clickedColumn = int((x - self.timeLineXOrigin) / self.offset)
            # print("clickedRow: " + str(self.clickedRow))
            self.rendered = False
            self.clickedRowSignal.emit()
            self.repaint()
            self.update()

    def mouseReleaseEvent(self, event):
        pass

    def mouseMoveEvent(self, event):
        # print(event.pos())
        pos = event.pos()
        x, y = pos.x(), pos.y()
        self.clickedColumn = -1
        if self.vismode and self.timeLineXOrigin < x < self.endOfTimeLine:
            self.clickedColumn = int((x - self.timeLineXOrigin) / self.offset)
            # print("clicked on frame %d",self.clickedColumn)
            self.clickedColumnSignal.emit()

    def switchToVisMode(self, vismode):
        """Visualize contacts directly in VMD by selecting a specific row."""
        self.vismode = vismode
        # self.labelView.vismode = vismode

    def paintEvent(self, event):

        qp = QPainter()
        qp.begin(self)

        # render pixmap to resolve performance issues
        if self.rendered:
            self.drawRenderedContact(qp)
        elif self.rendered is False and self.contacts:
            self.renderContact(False)
            self.rendered = True

        self.setMinimumSize(QSize(self.sizeX, self.sizeY))

        qp.end()

    def renderContact(self, generator):
        """Render the contact with the defined colors."""
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

        # probably included in next version...
        # if self.vismode:
        #     textoffset += 15
        #     start_text += 15
        #     startx += 15
        #     orig_startx += 15

        # print(orig_startx)
        self.timeLineXOrigin = orig_startx
        self.rowh = rowheight

        # self.sizeX = (len(self.contacts[0].scoreArray) + startx) * offset
        # self.sizeY = len(self.contacts) * rowheight
        if self.rangeFilterActive:
            self.sizeX = startx + (len(self.contacts[0].scoreArray) + merge * 2) * offset / merge
        else:
            self.sizeX = startx + (len(self.contacts[0].scoreArray[self.range[0]:self.range[1]]) + merge * 2) \
                                  * offset / merge

        # add one row for frame numbers
        self.sizeY = (len(self.contacts)+1) * rowheight

        self.pixmap = QPixmap(QSize(self.sizeX, self.sizeY))
        p = QPainter()

        if generator:
            p.begin(generator)
        else:
            p.begin(self.pixmap)

        p.fillRect(0, 0, self.sizeX, self.sizeY, whiteColor)

        row = 0
        rownumber = 0
        # print("merge value", merge)
        for c in self.contacts:
            self.alphaFactor = 50
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
            if rownumber == 0:
                # show the frame numbers on top
                p.setFont(QFont('Arial', 8))
                p.drawText(start_text, row + textoffset + 2.0, "Frame:")

                off = 0
                if self.range[0] != 0:
                    off = 1
                for l in range(self.range[0] + off, self.range[1] + 1, 10)[off:]:
                    if l == 0:
                        continue
                    # print(l)
                    # TODO: sometimes errors occur!
                    p.drawText(startx + (l - 1 - self.range[0]) * offset, row + textoffset + 2.0, str(l * merge))
                self.labelView.move(0, rowheight)
                row += rowheight
            while i < len(rangedScores):
                p.setPen(blackColor)
                merged_score = 0
                for j in range(merge):
                    if (i + j) >= len(rangedScores):
                        break
                    x = rangedScores[i + j]
                    merged_score += x
                merged_score /= merge
                alpha = merged_score * self.alphaFactor
                if alpha > 255:
                    alpha = 255
                if self.colorScheme == ColorScheme.bbsc:
                    # pass
                    p.setBrush(QColor(bbScColor[0], bbScColor[1], bbScColor[2], alpha))
                elif self.colorScheme == ColorScheme.custom:
                    color = QColor(self.customColor)
                    color.setAlpha(alpha)
                    p.setBrush(color)

                if rownumber == self.clickedRow:
                    p.setPen(QColor(250, 50, 50))
                else:
                    p.setPen(QColor(0, 0, 0))
                p.drawRect(startx, row, offset, 20)
                startx += offset
                i += merge
            self.offset = offset
            self.endOfTimeLine = startx
            startx = orig_startx
            row += rowheight
            rownumber += 1

        if generator:
            row = rowheight
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
        self.globalClickedRow = self.clickedRow
        self.clickedRow = -1

    def drawRenderedContact(self, qp):
        """Draws the rendered contact to the canvas."""
        qp.drawPixmap(0, 0, self.sizeX, self.sizeY, self.pixmap)
