import numpy as np

from .multi_trajectory import run_load_parallel
from .Biochemistry import (AccumulatedContact, AtomContact,
                           AccumulationMapIndex, AtomType, HydrogenBond,
                           AtomHBondType, TempContactAccumulate,
                           HydrogenBondAtoms)
from .ContactTrajectories import (AtomicContactTrajectory,
                                  AccumulatedContactTrajectory)
from .KeyManager import KeyManager

class ContactAnalyzer:

    def __init__(self, atomicTrajectory, map1, map2):
        self.resname_array = atomicTrajectory.resname_array
        self.resid_array = atomicTrajectory.resid_array
        self.name_array = atomicTrajectory.name_array
        self.segids = atomicTrajectory.segids
        self.backbone = atomicTrajectory.backbone
        self.atomicContacts = atomicTrajectory.contacts
        self.map1, self.map2 = map1, map2

    def accumulateContacts(self):
        numberOfFrames = len(self.atomicContacts)
        contactScores = np.zeros([0, numberOfFrames])
        hbonds = np.zeros([0, numberOfFrames])
        bbScores1 = np.array([])
        bbScores2 = np.array([])
        scScores1 = np.array([])
        scScores2 = np.array([])
        keys = np.array([], dtype=object)

        keyManager = KeyManager(self.name_array,
                                self.resid_array,
                                self.resname_array,
                                self.segids)

        for frame_id, frame_data in enumerate(self.atomicContacts):
            for contact in frame_data:
                key1, key2 = keyManager.makeKeyArraysFromMaps(self.map1,
                                                              self.map2,
                                                              contact)
                key = KeyManager.makeKeyFromKeyArrays(key1, key2)
                searchResult = np.where(keys == key)[0]

                contactIndex = -1
                if len(searchResult) != 0:
                    contactIndex = searchResult[0]
                else:
                    keys = np.append(keys, key)
                    bbScores1 = np.append(bbScores1, 0)
                    bbScores2 = np.append(bbScores2, 0)
                    scScores1 = np.append(scScores1, 0)
                    scScores2 = np.append(scScores2, 0)
                    contactScores = np.vstack((contactScores, np.zeros(numberOfFrames)))
                    hbonds = np.vstack((hbonds, np.zeros(numberOfFrames)))

                currentWeight = contact.weight
                contactScores[contactIndex, frame_id] += currentWeight
                hbonds[contactIndex, frame_id] += len(contact.hbondinfo)
                if contact.idx1 in self.backbone:
                    bbScores1[contactIndex] += currentWeight
                else:
                    scScores1[contactIndex] += currentWeight
                if contact.idx2 in self.backbone:
                    bbScores2[contactIndex] += currentWeight
                else:
                    scScores2[contactIndex] += currentWeight



        # print(np.count_nonzero(hbonds, axis=1)/numberOfFrames)
        ct = AccumulatedContactTrajectory(keys=keys,
                               contactScores=contactScores,
                               bbScores=[bbScores1, bbScores2],
                               scScores=[scScores1, scScores2],
                               hbonds=hbonds)
        return ct


    @staticmethod
    def runFrameScan(topology, trajectory, trajectoryScanParameters, nproc, use_pmda=True):
        if use_pmda:
            # TODO: refactor
            import MDAnalysis as mda
            from .ParallelAnalysis import AtomicContacts
            u = mda.Universe(topology, trajectory)
            print("starting analysis with pmda")
            ac = AtomicContacts(u,
                                trajectoryScanParameters.sel1text,
                                trajectoryScanParameters.sel2text,
                                [trajectoryScanParameters.cutoff,
                                trajectoryScanParameters.hbondcutoff,
                                trajectoryScanParameters.hbondcutangle])
            ac.run(n_blocks=nproc)
            return AtomicContactTrajectory(*ac.result)
        else:
            print("starting analysis with multiprocessing")
            return AtomicContactTrajectory(*run_load_parallel(nproc, topology,
                                                             trajectory,
                                                             trajectoryScanParameters.cutoff,
                                                             trajectoryScanParameters.hbondcutoff,
                                                             trajectoryScanParameters.hbondcutangle,
                                                             trajectoryScanParameters.sel1text,
                                                             trajectoryScanParameters.sel2text))
