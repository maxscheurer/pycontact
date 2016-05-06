from PyQt5.QtWidgets import QApplication
from MainWindow import *
def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()