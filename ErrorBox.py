from PyQt5.QtWidgets import QMessageBox

class ErrorBox(QMessageBox):
    """docstring for ErrorBox"""
    def __init__(self, msg):
        super(ErrorBox, self).__init__()
        self.msg = msg
        self.setIcon(QMessageBox.Warning)
        self.setText(self.msg)
        self.setWindowTitle("Error")
        self.setStandardButtons(QMessageBox.Ok)
