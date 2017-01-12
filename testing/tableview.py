import operator
from PyQt5.QtCore import QAbstractTableModel,Qt,pyqtSlot
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
class MyWindow(QWidget):
    def __init__(self, data_list, header, *args):
        QWidget.__init__(self, *args)
        # setGeometry(x_pos, y_pos, width, height)
        self.setGeometry(300, 200, 570, 450)
        self.setWindowTitle("Click on column title to sort")
        table_model = MyTableModel(self, data_list, header)
        table_view = QTableView()
        table_view.setModel(table_model)
        # set font
        font = QFont("Courier New", 14)
        table_view.setFont(font)
        # set column width to fit contents (set font first!)
        table_view.resizeColumnsToContents()
        # enable sorting
        table_view.setSortingEnabled(True)
        layout = QVBoxLayout(self)
        layout.addWidget(table_view)
        self.setLayout(layout)
class MyTableModel(QAbstractTableModel):
    def __init__(self, parent, mylist, header, *args):
        QAbstractTableModel.__init__(self, parent, *args)
        self.mylist = mylist
        self.header = header
    def rowCount(self, parent):
        return len(self.mylist)
    def columnCount(self, parent):
        return len(self.mylist[0])
    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != Qt.DisplayRole:
            return None
        return self.mylist[index.row()][index.column()]
    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.header[col]
        return None

    @pyqtSlot()
    def sort(self, col, order):
        """sort table by given column number col"""
        # self.emit(SIGNAL("layoutAboutToBeChanged()"))
        self.layoutAboutToBeChanged()
        self.mylist = sorted(self.mylist,
            key=operator.itemgetter(col))
        if order == Qt.DescendingOrder:
            self.mylist.reverse()
        # self.emit(SIGNAL("layoutChanged()"))
        self.layoutChanged()
# the solvent data ...
header = ['Solvent Name', ' BP (deg C)', ' MP (deg C)', ' Density (g/ml)']
# use numbers for numeric data to sort properly
data_list = [
('ACETIC ACID', 117.9, 16.7, 1.049),
('ACETIC ANHYDRIDE', 140.1, -73.1, 1.087),
('ACETONE', 56.3, -94.7, 0.791),
('ACETONITRILE', 81.6, -43.8, 0.786),
('ANISOLE', 154.2, -37.0, 0.995),
('BENZYL ALCOHOL', 205.4, -15.3, 1.045),
('BENZYL BENZOATE', 323.5, 19.4, 1.112),
('BUTYL ALCOHOL NORMAL', 117.7, -88.6, 0.81),
('BUTYL ALCOHOL SEC', 99.6, -114.7, 0.805),
('BUTYL ALCOHOL TERTIARY', 82.2, 25.5, 0.786),
('CHLOROBENZENE', 131.7, -45.6, 1.111),
('CYCLOHEXANE', 80.7, 6.6, 0.779),
('CYCLOHEXANOL', 161.1, 25.1, 0.971),
('CYCLOHEXANONE', 155.2, -47.0, 0.947),
('DICHLOROETHANE 1 2', 83.5, -35.7, 1.246),
('DICHLOROMETHANE', 39.8, -95.1, 1.325),
('DIETHYL ETHER', 34.5, -116.2, 0.715),
('DIMETHYLACETAMIDE', 166.1, -20.0, 0.937),
('DIMETHYLFORMAMIDE', 153.3, -60.4, 0.944),
('DIMETHYLSULFOXIDE', 189.4, 18.5, 1.102),
('DIOXANE 1 4', 101.3, 11.8, 1.034),
('DIPHENYL ETHER', 258.3, 26.9, 1.066),
('ETHYL ACETATE', 77.1, -83.9, 0.902),
('ETHYL ALCOHOL', 78.3, -114.1, 0.789),
('ETHYL DIGLYME', 188.2, -45.0, 0.906),
('ETHYLENE CARBONATE', 248.3, 36.4, 1.321),
('ETHYLENE GLYCOL', 197.3, -13.2, 1.114),
('FORMIC ACID', 100.6, 8.3, 1.22),
('HEPTANE', 98.4, -90.6, 0.684),
('HEXAMETHYL PHOSPHORAMIDE', 233.2, 7.2, 1.027),
('HEXANE', 68.7, -95.3, 0.659),
('ISO OCTANE', 99.2, -107.4, 0.692),
('ISOPROPYL ACETATE', 88.6, -73.4, 0.872),
('ISOPROPYL ALCOHOL', 82.3, -88.0, 0.785),
('METHYL ALCOHOL', 64.7, -97.7, 0.791),
('METHYL ETHYLKETONE', 79.6, -86.7, 0.805),
('METHYL ISOBUTYL KETONE', 116.5, -84.0, 0.798),
('METHYL T-BUTYL ETHER', 55.5, -10.0, 0.74),
('METHYLPYRROLIDINONE N', 203.2, -23.5, 1.027),
('MORPHOLINE', 128.9, -3.1, 1.0),
('NITROBENZENE', 210.8, 5.7, 1.208),
('NITROMETHANE', 101.2, -28.5, 1.131),
('PENTANE', 36.1, ' -129.7', 0.626),
('PHENOL', 181.8, 40.9, 1.066),
('PROPANENITRILE', 97.1, -92.8, 0.782),
('PROPIONIC ACID', 141.1, -20.7, 0.993),
('PROPIONITRILE', 97.4, -92.8, 0.782),
('PROPYLENE GLYCOL', 187.6, -60.1, 1.04),
('PYRIDINE', 115.4, -41.6, 0.978),
('SULFOLANE', 287.3, 28.5, 1.262),
('TETRAHYDROFURAN', 66.2, -108.5, 0.887),
('TOLUENE', 110.6, -94.9, 0.867),
('TRIETHYL PHOSPHATE', 215.4, -56.4, 1.072),
('TRIETHYLAMINE', 89.5, -114.7, 0.726),
('TRIFLUOROACETIC ACID', 71.8, -15.3, 1.489),
('WATER', 100.0, 0.0, 1.0),
('XYLENES', 139.1, -47.8, 0.86)
]
app = QApplication([])
win = MyWindow(data_list, header)
win.show()
app.exec_()
