""" Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
"""

import warnings
import sys

from PyQt5.QtWidgets import QApplication

from MainWindow import MainWindow

warnings.filterwarnings("ignore")


def main():
    """Pycontact main function."""
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()


if __name__ == '__main__':
    main()
