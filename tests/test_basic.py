from unittest import TestCase
import sys
from os import path
from PyContact.core.ContactAnalyzer import *
from PyContact.exampleData.datafiles import DCD, PSF, TPR, XTC
from PyContact.core.Biochemistry import AtomContact
from PyContact.core.multi_trajectory import weight_function
from PyContact.cy_modules.cy_gridsearch import cy_find_contacts

import MDAnalysis as mda
import multiprocessing

import numpy as np
import pandas as pd

from MDAnalysis.lib.distances import capped_distance

multiprocessing.log_to_stderr()

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))


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

    def test_distances(self):
        u = mda.Universe(self.psffile, self.dcdfile)
        sel1 = u.select_atoms("segid RN11")
        sel2 = u.select_atoms("segid UBQ")

        for ts in u.trajectory:
            pairs, distances = capped_distance(
                sel1.positions, sel2.positions,
                max_cutoff=5.0, min_cutoff=None,
                box=None, method=None, return_distances=True
            )
            df = pd.DataFrame(data=pairs, columns=['atom1', 'atom2'])
            df['dist'] = distances

            # compute pairs and distances using C code
            natoms1 = len(sel1.atoms)
            natoms2 = len(sel2.atoms)
            pos1 = np.array(np.reshape(sel1.positions, (1, natoms1 * 3)), dtype=np.float64)
            pos2 = np.array(np.reshape(sel2.positions, (1, natoms2 * 3)), dtype=np.float64)
            xyz1 = np.array(pos1, dtype=np.float32)
            xyz2 = np.array(pos2, dtype=np.float32)
            nbList1 = cy_find_contacts(xyz1, natoms1, xyz2, natoms2, 5.0)

            idx1 = 0
            pairs_list = []
            distances_list = []
            for atom1sNeighbors in nbList1:
                for idx2 in atom1sNeighbors:
                    dvec = pos1[0][3*idx1:3*idx1+3] - pos2[0][3*idx2:3*idx2+3]
                    distance = np.sqrt(dvec.dot(dvec))
                    pairs_list.append(np.array([idx1, idx2], dtype=np.int64))
                    distances_list.append(distance)
                idx1 += 1
            df_old = pd.DataFrame(data=np.array(pairs_list), columns=['atom1', 'atom2'])
            df_old['dist'] = np.array(distances_list)

            df.sort_values(by=['atom1', 'atom2'], inplace=True)
            df_old.sort_values(by=['atom1', 'atom2'], inplace=True)

            # compare results
            np.testing.assert_allclose(df['dist'], df_old['dist'], atol=1e-5)
            np.testing.assert_equal(df['atom1'].values, df_old['atom1'].values)
            np.testing.assert_equal(df['atom2'].values, df_old['atom2'].values)

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
