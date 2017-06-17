from unittest import TestCase
import sys
from os import path
from PyContact.core.ContactAnalyzer import *
from PyContact.exampleData.datafiles import DCD, PSF, TPR, XTC
import MDAnalysis as mda
import multiprocessing
multiprocessing.log_to_stderr()

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
# from PyContact.core.multi_trajectory import run_load_parallel


class PsfDcdReadingTest(TestCase):
    def setUp(self):
        self.dcdfile = DCD
        self.psffile = PSF
        self.tpr = TPR
        self.xtc = XTC

    def tearDown(self):
        del self.dcdfile
        del self.psffile

    def test_import_dcd_file(self):
        mda.Universe(self.psffile, self.dcdfile)

    def test_import_xtc_file(self):
        # seg_0_Protein_chain_U
        # seg_1_Protein_chain_R
        mda.Universe(self.tpr, self.xtc)

    def test_singleCore_analysis(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
        analyzer.runFrameScan(1)
        self.assertEqual(len(analyzer.contactResults), 50)
        map1 = [0, 0, 1, 1, 0]
        map2 = [0, 0, 1, 1, 0]
        analyzer.runContactAnalysis(map1, map2, 1)
        self.assertEqual(len(analyzer.finalAccumulatedContacts), 148)
        hbond_sum = 0
        for c in analyzer.finalAccumulatedContacts:
            hbond_sum += c.hbond_percentage()
        self.assertEqual(hbond_sum, 676.0)

    def test_trackMolecule(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11 and resid 85", "segid UBQ")
        analyzer.runFrameScan(1)
        analyzer.runMoleculeTracking(1, [0, 0, 1, 1, 0])

    def test_selfInteraction_analysis(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11", "self")
        analyzer.runFrameScan(1)
        self.assertEqual(len(analyzer.contactResults), 50)
        map1 = [0, 0, 1, 1, 0]
        map2 = [0, 0, 1, 1, 0]
        analyzer.runContactAnalysis(map1, map2, 1)
        # self.assertEqual(len(analyzer.finalAccumulatedContacts), 148)
        # hbond_sum = 0
        # for c in analyzer.finalAccumulatedContacts:
            # hbond_sum += c.hbond_percentage()
        # self.assertEqual(hbond_sum, 676.0)

    def test_zero_atomselection(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid A", "resid 100")
        try:
            analyzer.runFrameScan(1)
        except:
            print("Error in atom selection caught.")

        try:
            analyzer.runFrameScan(4)
        except:
            print("Error in atom selection (multicore) caught.")

    def test_selfInteraction_analysis_parallel(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11", "self")
        analyzer.runFrameScan(2)
        self.assertEqual(len(analyzer.contactResults), 50)
        map1 = [0, 0, 1, 1, 0]
        map2 = [0, 0, 1, 1, 0]
        analyzer.runContactAnalysis(map1, map2, 1)

    def test_multiCore_analysis(self):
        analyzer = Analyzer(self.psffile, self.dcdfile, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
        analyzer.runFrameScan(2)
        self.assertEqual(len(analyzer.contactResults), 50)
        map1 = [0, 0, 1, 1, 0]
        map2 = [0, 0, 1, 1, 0]
        analyzer.runContactAnalysis(map1, map2, 2)
        self.assertEqual(len(analyzer.finalAccumulatedContacts), 148)
        hbond_sum = 0
        for c in analyzer.finalAccumulatedContacts:
            hbond_sum += c.hbond_percentage()
        self.assertEqual(hbond_sum, 676.0)

    def test_around_selection_patch(self):
        univ = mda.Universe(self.psffile, self.dcdfile)
        aroundText = "segid UBQ and around 5 segid RN11"
        sel = univ.select_atoms(aroundText)
        self.assertEqual(len(sel), 261)
