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
        mda.Universe(self.psffile,self.dcdfile)

    def test_simple_analysis(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
        analyzer.runFrameScan()
        self.assertEqual(len(analyzer.contactResults), 50)
        map1=[0,0,0,1,1,0]
        map2=[0,0,0,1,1,0]
        analyzer.runContactAnalysis(map1, map2)
        self.assertEqual(len(analyzer.finalAccumulatedContacts), 148)
