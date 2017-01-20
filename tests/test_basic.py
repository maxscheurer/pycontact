from unittest import TestCase
import sys
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
from ..ContactAnalyzer import * 
import MDAnalysis as mda

class PsfDcdReadingTest(TestCase):
    def setUp(self):
        self.dcdfile = "../short.dcd"
        self.psffile = "../rpn11_ubq_interface-ionized.psf"

    def tearDown(self):
        del self.dcdfile
        del self.psffile

    def test_import_dcd_file(self):
        self.u = mda.Universe(self.psffile,self.dcdfile)
