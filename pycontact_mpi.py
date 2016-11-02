'''
    Authors: Maximilian Scheurer, Peter Rodenkirch
    Date created: May 2016
    Python Version: 2.7
    Version: 0.1a
    Status: Development
'''
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size


if rank == 0:
    from MainWindow import *
    def main():
        app = QApplication(sys.argv)
        window = MainWindow()
        window.show()
        app.exec_()

    if __name__ == '__main__':
        main()
