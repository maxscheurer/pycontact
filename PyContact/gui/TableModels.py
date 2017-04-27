from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QVariant
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTableView, QWidget, QVBoxLayout
import pandas as pd

class Widget(QWidget):
    """
    A simple test widget to contain and own the model and table.
    """
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        box = QVBoxLayout()
        self.setLayout(box)
        cdf = self.get_data_frame()
        self._tm = ResidueTableModel(self)
        self._tm.update(cdf)
        self._tv = TableView(self)
        self._tv.setModel(self._tm)
        box.addWidget(self._tv)

    def get_data_frame(self):
        df = pd.DataFrame({'Name':['a','b','c','d'],
        'First':[2.3,5.4,3.1,7.7], 'Last':[23.4,11.2,65.3,88.8], 'Class':[1,1,2,1], 'Valid':[True, True, True, False]})
        return df

class TableView(QTableView):
    """
    A simple table to demonstrate the QComboBox delegate.
    """
    def __init__(self, *args, **kwargs):
        QTableView.__init__(self, *args, **kwargs)

class ResidueTableModel(QAbstractTableModel):
    def __init__(self, parent=None, *args):
        super(ResidueTableModel, self).__init__()
        self.datatable = None

    def update(self, dataIn):
        self.datatable = dataIn
        print('Datatable : {0}'.format(self.datatable))

    def rowCount(self, parent=QModelIndex()):
        return len(self.datatable.index)

    def columnCount(self, parent=QModelIndex()):
        return len(self.datatable.columns.values)

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            i = index.row()
            j = index.column()
            return '{0}'.format(self.datatable.iget_value(i, j))
        else:
            return QVariant()

    def flags(self, index):
        return Qt.ItemIsEnabled
