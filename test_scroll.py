from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import (QWidget, QApplication)
from PyQt5.QtWidgets import (QApplication, QWidget, QDesktopWidget, 
    QLabel, QCheckBox, QPushButton, QMainWindow, QMenuBar, QComboBox,
    QLineEdit, QTextEdit, QGridLayout, QFileDialog, QAction, qApp, QHBoxLayout, QGroupBox, QFormLayout, QScrollArea,
    QVBoxLayout)
class Window(QWidget):
    def __init__(self, val):
        QWidget.__init__(self)
        mygroupbox = QGroupBox('this is my groupbox')
        myform = QFormLayout()
        labellist = []
        combolist = []
        for i in range(val):
            labellist.append(QLabel('mylabel'))
            combolist.append(QComboBox())
            myform.addRow(labellist[i],combolist[i])
        mygroupbox.setLayout(myform)
        scroll = QScrollArea()
        scroll.setWidget(mygroupbox)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QVBoxLayout(self)
        layout.addWidget(scroll)

if __name__ == '__main__':

    import sys
    app = QApplication(sys.argv)
    window = Window(25)
    window.setGeometry(500, 300, 300, 400)
    window.show()
    sys.exit(app.exec_())
