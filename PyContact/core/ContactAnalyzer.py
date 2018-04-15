import numpy as np

from .multi_trajectory import run_load_parallel
from .Biochemistry import (AccumulatedContact, AtomContact,
                           AccumulationMapIndex, AtomType, HydrogenBond,
                           AtomHBondType, TempContactAccumulate,
                           HydrogenBondAtoms)
from .ContactTrajectories import (AtomicContactTrajectory,
                                  AccumulatedContactTrajectory)

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
        keys = np.array([])

        for frame_id, frame_data in enumerate(self.atomicContacts):
            for contact in frame_data:
                key1, key2 = self.makeKeyArraysFromMaps(self.map1,
                                                        self.map2,
                                                        contact)
                key = self.makeKeyFromKeyArrays(key1, key2)
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
        np.save("keys.np", keys)
        np.save("contactScores.np", contactScores)
        ct = AccumulatedContactTrajectory(keys=keys,
                               contactScores=contactScores,
                               bbScores=[bbScores1, bbScores2],
                               scScores=[scScores1, scScores2],
                               hbonds=hbonds)
        return ct

    def makeKeyArraysFromMaps(self, map1, map2, contact):
        """Creates key Arrays from the chosen accumulation maps.

            maps contain information whether to consider an atom's property for contact accumulation
            map1 and map2 contain 5 boolean values each, cf. AccumulationMapIndex
            for a given contact, the corresponding value to a property is written to keys1 and keys2, respectively
            example input:
            map1 = [0,0,1,1,0]
            map2 = [0,0,1,1,0], meaning that residue and resname should be used for contact accumulation
            contact: idx1,idx2
            results: (example!)
            keys1=["none","none","none","14", "VAL", "none"]
            keys2=["none","none","none","22", "ILE, "none"]

        """
        idx1 = contact.idx1
        idx2 = contact.idx2
        counter = 0
        keys1 = []
        for val in map1:
            if val == 1:
                if counter == AccumulationMapIndex.index:
                    keys1.append(idx1)
                elif counter == AccumulationMapIndex.name:
                    keys1.append(self.name_array[idx1])
                elif counter == AccumulationMapIndex.resid:
                    keys1.append(self.resid_array[idx1])
                elif counter == AccumulationMapIndex.resname:
                    keys1.append(self.resname_array[idx1])
                elif counter == AccumulationMapIndex.segid:
                    keys1.append(self.segids[idx1])
            else:
                keys1.append("none")
            counter += 1
        counter = 0
        keys2 = []
        for val in map2:
            if val == 1:
                if counter == AccumulationMapIndex.index:
                    keys2.append(idx2)
                elif counter == AccumulationMapIndex.name:
                    keys2.append(self.name_array[idx2])
                elif counter == AccumulationMapIndex.resid:
                    keys2.append(self.resid_array[idx2])
                elif counter == AccumulationMapIndex.resname:
                    keys2.append(self.resname_array[idx2])
                elif counter == AccumulationMapIndex.segid:
                    keys2.append(self.segids[idx2])
            else:
                keys2.append("none")
            counter += 1
        return [keys1, keys2]

    @staticmethod
    def makeKeyFromKeyArrays(key1, key2):
        """Returns a human readable key from two key arrays.
            example:
            keys1=["none","none","14", "VAL", "none"]
            keys2=["none","none","22", "ILE", "none"]
            returns a human readable key with the mapping identifiers in AccumulationMapIndex
            in the given example data:
            key="r.14rn.VAL-r.22rn.ILE"
            key is used to accumulated contacts in a dictionary (= a contact's unique identifier)
        """
        key = ""
        itemcounter = 0
        for item in key1:
            if item != "none":
                key += AccumulationMapIndex.mapping[itemcounter] + str(item)
            itemcounter += 1
        key += "-"
        itemcounter = 0
        for item in key2:
            if item != "none":
                key += AccumulationMapIndex.mapping[itemcounter] + str(item)
            itemcounter += 1
        return key


    @staticmethod
    def runFrameScan(topology, trajectory, trajectoryScanParameters, nproc):
        return AtomicContactTrajectory(*run_load_parallel(nproc, topology,
                                                         trajectory, trajectoryScanParameters.cutoff,
                                                         trajectoryScanParameters.hbondcutoff,
                                                         trajectoryScanParameters.hbondcutangle,
                                                         trajectoryScanParameters.sel1text,
                                                         trajectoryScanParameters.sel2text))
