from PyQt5.QtWidgets import QAbstractButton
from PyQt5.QtGui import QPixmap, QImage, QPainter


class PicButton(QAbstractButton):
    def __init__(self, pixmap, parent=None):
        super(PicButton, self).__init__(parent)
        self.pixmap = pixmap
        self.setMaximumWidth(30)
        self.setMaximumHeight(30)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.drawPixmap(event.rect(), self.pixmap)

    def sizeHint(self):
        return self.pixmap.size()


class HelpButton(PicButton):
    def __init__(self, parent=None):
        helpImage = QImage("/home/max/help.png")
        helpPixmap = QPixmap.fromImage(helpImage)
        super(HelpButton, self).__init__(helpPixmap, parent)
