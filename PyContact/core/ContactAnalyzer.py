from __future__ import print_function
import itertools
import re
import os
import time
import multiprocessing
from multi_accumulation import *
from multi_trajectory import run_load_parallel
from LogPool import *
from copy import deepcopy
import operator

import pickle

import MDAnalysis
from MDAnalysis.analysis import distances
import numpy as np
from PyQt5.QtCore import pyqtSlot, pyqtSignal, QObject

# TODO: fix aroundPatch with gridsearch in C code using cython
from .aroundPatch import AroundSelection
from .Biochemistry import (AccumulatedContact, AtomContact, AccumulationMapIndex, AtomType, HydrogenBond, AtomHBondType, TempContactAccumulate, HydrogenBondAtoms)
from ..cy_modules.cy_gridsearch import cy_find_contacts

MDAnalysis.core.flags['use_periodic_selections'] = False
MDAnalysis.core.flags['use_KDTree_routines'] = True


class Analyzer(QObject):
    """Performs a contact search and analyzes the results."""
    frameUpdate = pyqtSignal(float)

    def __init__(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        super(Analyzer, self).__init__()
        self.psf = psf
        self.dcd = dcd
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.sel1text = sel1text
        self.sel2text = sel2text
        self.lastMap1 = []
        self.lastMap2 = []
        self.contactResults = []
        self.contactResults = []
        self.resname_array = []
        self.resid_array = []
        self.name_array = []
        self.segids = []
        self.backbone = []
        self.finalAccumulatedContacts = []
        self.bonds = []
        self.totalFrameNumber = 0
        self.currentFrameNumber = 0
        self.analysis_state = False
        self.totalFramesToProcess = 1

    def runFrameScan(self, nproc):
        """Performs a contact search using nproc threads."""
        try:
            if nproc == 1:
                self.contactResults = self.analyze_psf_dcd_grid(self.psf, self.dcd, self.cutoff, self.hbondcutoff,
                                                       self.hbondcutangle, self.sel1text, self.sel2text)
                # old: original version but slower
                # self.contactResults = self.analyze_psf_dcd(self.psf, self.dcd, self.cutoff, self.hbondcutoff,
                #                                        self.hbondcutangle, self.sel1text, self.sel2text)
            else:
                self.contactResults, self.resname_array, self.resid_array, self.name_array, self.segids, \
                            self.backbone = run_load_parallel(nproc, self.psf, self.dcd, self.cutoff, self.hbondcutoff,
                                                              self.hbondcutangle, self.sel1text, self.sel2text)
        except:
            raise Exception

    def runContactAnalysis(self, map1, map2, nproc):
        """Performs a contadt analysis using nproc threads."""
        if nproc == 1:
            self.finalAccumulatedContacts = self.analyze_contactResultsWithMaps(self.contactResults, map1, map2)
        else:
            self.finalAccumulatedContacts = self.analyze_contactResultsWithMaps_Parallel(map1, map2, nproc)

        self.lastMap1 = map1
        self.lastMap2 = map2
        return deepcopy(self.finalAccumulatedContacts)

    def runMoleculeTracking(self, selIndex, map):
        print("Running tracking")
        return self.analyze_trackMolecule(self.contactResults, selIndex, map)

    def setTrajectoryData(
            self, resname_array, resid_array, name_array,
            segids, backbone,
            sel1text, sel2text):
        self.resname_array = resname_array
        self.resid_array = resid_array
        self.name_array = name_array
        self.segids = segids
        self.backbone = backbone
        self.sel1text = sel1text
        self.sel2text = sel2text

    def getTrajectoryData(self):
        return [self.resname_array, self.resid_array, self.name_array,
                self.segids, self.backbone, self.sel1text,
                self.sel2text]

    def getFilePaths(self):
        return self.psf, self.dcd

    # find a string in s between the strings first and last
    @staticmethod
    def find_between(s, first, last):
        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

    @staticmethod
    def weight_function(value):
        """weight function to score contact distances"""
        return 1.0 / (1.0 + np.exp(5.0 * (value - 4.0)))

    def makeKeyArraysFromMaps(self, map1, map2, contact):
        """Creates key Arrays from the chosen accumulation maps.

            maps contain information wether to consider an atom's field for contact accumulation
            map1 and map2 contain six boolean values each, cf. AccumulationMapIndex
            for a given contact, the corresponding value to a field is written to keys1 and keys2, respectively
            example input:
            map1 = [0,0,0,1,1,0]
            map2 = [0,0,0,1,1,0], meaning that residue and resname should be used for contact accumulation
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

    def makeKeyArraysFromKey(self, key):
        """Converts a key to two key arrays.
            "inverse" function of makeKeyFromKeyArrays
        """
        keystring1, keystring2 = key.split("-")
        mapping = AccumulationMapIndex.mapping
        maximal = len(mapping)
        key1 = []
        for i in range(0, maximal):
            current = mapping[i]
            if current not in keystring1:
                key1.append("none")
                continue
            if i == (maximal - 1):
                key1.append(keystring1[keystring1.index(current) + len(current):])
                break
            nextCurrent = mapping[i + 1]
            if nextCurrent not in keystring1:
                nxt = ""
                for k in mapping[i + 1:]:
                    if k in keystring1[keystring1.index(current) + len(current):]:
                        nxt = k
                        break
                if nxt != "":
                    key1.append(keystring1[keystring1.index(current) + len(current):keystring1.index(nxt)])
                else:
                    key1.append(keystring1[keystring1.index(current) + len(current):])
                continue
            else:
                currentValue = self.find_between(keystring1, current, nextCurrent)
                if currentValue == "":
                    key1.append("none")
                else:
                    key1.append(currentValue)

        key2 = []
        for i in range(0, maximal):
            current = mapping[i]
            if current not in keystring2:
                key2.append("none")
                continue
            if i == (maximal - 1):
                key2.append(keystring2[keystring2.index(current) + len(current):])
                break
            nextCurrent = mapping[i + 1]
            if nextCurrent not in keystring2:
                nxt = ""
                for k in mapping[i + 1:]:
                    if k in keystring2[keystring2.index(current) + len(current):]:
                        nxt = k
                        break
                if nxt != "":
                    key2.append(keystring2[keystring2.index(current) + len(current):keystring2.index(nxt)])
                else:
                    key2.append(keystring2[keystring2.index(current) + len(current):])
                continue
            else:
                currentValue = self.find_between(keystring2, current, nextCurrent)
                if currentValue == "":
                    key2.append("none")
                else:
                    key2.append(currentValue)
        return [key1, key2]

    @staticmethod
    def make_single_title(key):
        """returns the title of the AccumulatedContact to be displayed in contact's label"""
        titleDict = {}
        counter = 0
        for item in key:
            titleDict[AccumulationMapIndex.mapping[counter]] = (item if item != "none" else "")
            counter += 1
        residueString = "%s%s" % (titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.resname]],
                                  str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.resid]]))
        atomIndexString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.index],
                                      str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.index]])) if
                           titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.index]] != "" else "")
        atomNameString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.name],
                                     str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.name]])) if
                          titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.name]] != "" else "")
        segnameString = ("%s %s" % (AccumulationMapIndex.mapping[AccumulationMapIndex.segid],
                                    str(titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.segid]])) if
                         titleDict[AccumulationMapIndex.mapping[AccumulationMapIndex.segid]] != "" else "")
        tempList = [residueString, atomIndexString, atomNameString, segnameString]
        finishedList = []
        for string in tempList:
            if string != "":
                finishedList.append(string)
        finishedString = ", ".join(finishedList)
        return finishedString

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


    # newer version, uses gridsearch to find contacts, much faster and less memory-intense
    def analyze_psf_dcd_grid(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        """Reading topology/trajectory and assessing hbonds"""

        # load psf and dcd file in memory
        u = MDAnalysis.Universe(psf, dcd)

# TODO: think about doing u.select_atoms(sel1text + " or " + sel2text)
        all_sel = u.select_atoms("all")
        # all_sel = u.select_atoms("%s or %s or name H.*" % (sel1text, sel2text))
        backbone_sel = u.select_atoms("backbone")
        self.resname_array = []
        self.resid_array = []
        self.name_array = []
        self.segids = []
        self.backbone = []
        for atom in all_sel.atoms:
            self.resname_array.append(atom.resname)
            self.resid_array.append(atom.resid)
            self.name_array.append(atom.name)
            self.bonds.append(atom.bonds)
            self.segids.append(atom.segid)
        for atom in backbone_sel:
            self.backbone.append(atom.index)

        selfInteraction = False

        if sel2text == "self":
            sel1 = u.select_atoms(sel1text)
            sel2 = u.select_atoms(sel1text)
            selfInteraction = True
        else:
            sel1 = u.select_atoms(sel1text)
            sel2 = u.select_atoms(sel2text)

        if (len(sel1.atoms) == 0 or len(sel2.atoms) == 0):
            raise Exception

        contactResults = []
        # loop over trajectory
        self.totalFrameNumber = len(u.trajectory)
        start = time.time()
        for ts in u.trajectory:
            # define selections according to sel1text and sel2text
            if "around" in sel1text:
                sel1 = u.select_atoms(sel1text)
            if "around" in sel2text:
                sel2 = u.select_atoms(sel2text)
            # write atomindices for each selection to list
            indices1 = []
            for at in sel1.atoms:
                indices1.append(at.index)
            indices2 = []
            for at in sel2.atoms:
                indices2.append(at.index)

            currentFrameContacts = []
            frame = ts.frame
            self.currentFrameNumber = ts.frame

            # pass positions and distance to cython
            natoms1 = len(sel1.atoms)
            natoms2 = len(sel2.atoms)
            pos1 = np.reshape(sel1.positions, (1, natoms1 * 3))
            pos2 = np.reshape(sel2.positions, (1, natoms2 * 3))
            xyz1 = np.array(pos1, dtype=np.float32)
            xyz2 = np.array(pos2, dtype=np.float32)
            # 2d array with index of atom1 being the index of the first dimension
            # individual lists contain atom2 indices
            nbList1 = cy_find_contacts(xyz1, natoms1, xyz2, natoms2, cutoff)

            # we only need the 1st list
            # nbList1 = res[:natoms1]
            # nbList2 = res[natoms1:]


            idx1 = 0
            for atom1sNeighbors in nbList1:
                for idx2 in atom1sNeighbors:
                    convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
                    convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
                    # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered,
                    # hydrogen bonds can still be detected!
                    if re.match("H(.*)", self.name_array[convindex1]) or re.match("H(.*)", self.name_array[convindex2]):
                        continue
                        # distance between atom1 and atom2
                    # check if residues are more than 4 apart, and in the same segment
                    if selfInteraction:
                        if (self.resid_array[convindex1] - self.resid_array[convindex2]) < 5 and self.segids[convindex1] == self.segids[convindex2]:
                            continue
                    # distance = distarray[idx1, idx2]
                    # weight = self.weight_function(distance)
                    distance = np.linalg.norm(pos1[0][3*idx1:3*idx1+3] - pos2[0][3*idx2:3*idx2+3])
                    # if (distance - distarray[idx1, idx2]) > 0.001:
                    #     print("Error in distance calculations!")
                    #     return
                    # # print(convindex1, convindex2, distance, distarray[idx1, idx2])
                    # if (distance > cutoff):
                    #     print("Distances must be smaller/equal cutoff!")
                    #     return
                    weight = self.weight_function(distance)

                    # HydrogenBondAlgorithm
                    hydrogenBonds = []
                    # FF independent hydrogen bonds
                    if (self.name_array[convindex1][0] in HydrogenBondAtoms.atoms and self.name_array[convindex2][0] in HydrogenBondAtoms.atoms):
                            # print("hbond? %s - %s" % (type_array[convindex1], type_array[convindex2]))
                            # search for hatom, check numbering in bond!!!!!!!!!!
                            b1 = self.bonds[convindex1]
                            b2 = self.bonds[convindex2]

                            # b1 = all_sel[convindex1].bonds
                            # b2 = all_sel[convindex2].bonds
                            # search for hydrogen atoms bound to atom 1
                            bondcount1 = 0
                            hydrogenAtomsBoundToAtom1 = []

                            # new code
                            for bnd in b1:
                                b = bnd.type
                                hydrogen = next((x for x in b if x.startswith("H")), 0)
                                # print(b)
                                if hydrogen != 0:
                                    # print("h bond to atom1")
                                    bondindices1 = b1.to_indices()[bondcount1]
                                    # print bondindices1
                                    # for j in bondindices1:
                                    #     print(self.type_array[j+1])
                                    hydrogenidx = next(
                                        (j for j in bondindices1 if self.name_array[j].startswith("H")), -1)
                                    if hydrogenidx != -1:
                                        # print(self.type_array[hydrogenidx])
                                        hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                                bondcount1 += 1
                            # search for hydrogen atoms bound to atom 2
                            bondcount2 = 0
                            hydrogenAtomsBoundToAtom2 = []
                            # print(b2)
                            for bnd2 in b2:
                                b = bnd2.type
                                hydrogen = next((x for x in b if x.startswith("H")), 0)
                                # print(b)
                                if hydrogen != 0:
                                    # print("h bond to atom2")
                                    bondindices2 = b2.to_indices()[bondcount2]
                                    hydrogenidx = next(
                                        (k for k in bondindices2 if self.name_array[k].startswith("H")), -1)
                                    if hydrogenidx != -1:
                                        # print(type_array[hydrogenidx])
                                        hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                                bondcount2 += 1
                            # check hbond criteria for hydrogen atoms bound to first atom
                            for global_hatom in hydrogenAtomsBoundToAtom1:
                                conv_hatom = indices1.index(global_hatom)
                                # print(typeHeavy)
                                #
                                # TODO: FF independent version
                                # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[conv_hatom, idx2] <= hbondcutoff):
                                # dist = distarray[conv_hatom, idx2]
                                # dist = np.linalg.norm(sel1.positions[conv_hatom] - sel2.positions[idx2])
                                dist = np.linalg.norm(pos1[0][3*conv_hatom:3*conv_hatom+3] - pos2[0][3*idx2:3*idx2+3])
                                if (dist <= hbondcutoff):
                                    donorPosition = sel1.positions[idx1]
                                    hydrogenPosition = sel1.positions[conv_hatom]
                                    acceptorPosition = sel2.positions[idx2]
                                    v1 = hydrogenPosition - acceptorPosition
                                    v2 = hydrogenPosition - donorPosition
                                    v1norm = np.linalg.norm(v1)
                                    v2norm = np.linalg.norm(v2)
                                    dot = np.dot(v1, v2)
                                    angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                                    # print(angle)
                                    if angle >= hbondcutangle:
                                        # print("new hbond")
                                        new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                                 hbondcutoff,
                                                                 hbondcutangle)
                                        hydrogenBonds.append(new_hbond)
                                        # print(str(convindex1) + " " + str(convindex2)
                                        # print("hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2))
                                        # print(angle)
                            for global_hatom in hydrogenAtomsBoundToAtom2:
                                conv_hatom = indices2.index(global_hatom)
                                # TODO: FF independent version
                                # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[idx1, conv_hatom] <= hbondcutoff):
                                # FIXME: WTF?
                                # if (distarray[conv_hatom, idx2] <= hbondcutoff):
                                # dist = distarray[idx1, conv_hatom]
                                # dist = np.linalg.norm(sel1.positions[idx1] - sel2.positions[conv_hatom])
                                dist = np.linalg.norm(pos1[0][3*idx1:3*idx1+3] - pos2[0][3*conv_hatom:3*conv_hatom+3])
                                if (dist <= hbondcutoff):
                                    donorPosition = sel2.positions[idx2]
                                    hydrogenPosition = sel2.positions[conv_hatom]
                                    acceptorPosition = sel1.positions[idx1]
                                    v1 = hydrogenPosition - acceptorPosition
                                    v2 = hydrogenPosition - donorPosition
                                    v1norm = np.linalg.norm(v1)
                                    v2norm = np.linalg.norm(v2)
                                    dot = np.dot(v1, v2)
                                    angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                                    if angle >= hbondcutangle:
                                        new_hbond = HydrogenBond(convindex2, convindex1, global_hatom, dist, angle,
                                                                 hbondcutoff,
                                                                 hbondcutangle)
                                        hydrogenBonds.append(new_hbond)
                                        # print str(convindex1) + " " + str(convindex2)
                                        # print "hbond found: %d,%d,%d"%(convindex2,global_hatom,convindex1)
                                        # print angle
                                        # finalize
                    newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),
                                                 int(convindex2),
                                                 hydrogenBonds)
                    currentFrameContacts.append(newAtomContact)
                idx1 += 1
            contactResults.append(currentFrameContacts)

        # pickle.dump(contactResults, open("single_results_experimental.dat", "w"))
        # for f in contactResults:
        #     print("experiment", len(f))
        stop = time.time()
        print("grid:",stop-start)
        return contactResults



    def analyze_psf_dcd(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        """Reading topology/trajectory and assessing hbonds"""

        # load psf and dcd file in memory
        u = MDAnalysis.Universe(psf, dcd)

# TODO: think about doing u.select_atoms(sel1text + " or " + sel2text)
        all_sel = u.select_atoms("all")
        # all_sel = u.select_atoms("%s or %s or name H.*" % (sel1text, sel2text))
        backbone_sel = u.select_atoms("backbone")
        self.resname_array = []
        self.resid_array = []
        self.name_array = []
        self.segids = []
        self.backbone = []
        for atom in all_sel.atoms:
            self.resname_array.append(atom.resname)
            self.resid_array.append(atom.resid)
            self.name_array.append(atom.name)
            self.bonds.append(atom.bonds)
            self.segids.append(atom.segid)
        for atom in backbone_sel:
            self.backbone.append(atom.index)


        # check if self-interaction is wanted
        selfInteraction = False

        if sel2text == "self":
            sel1 = u.select_atoms(sel1text)
            sel2 = u.select_atoms(sel1text)
            selfInteraction = True
        else:
            sel1 = u.select_atoms(sel1text)
            sel2 = u.select_atoms(sel2text)

        if (len(sel1.atoms) == 0 or len(sel2.atoms) == 0):
            raise Exception

        contactResults = []
        # loop over trajectory
        self.totalFrameNumber = len(u.trajectory)
        start = time.time()
        for ts in u.trajectory:
            # define selections according to sel1text and sel2text
            if "around" in sel1text:
                sel1 = u.select_atoms(sel1text)
            if "around" in sel2text:
                sel2 = u.select_atoms(sel2text)
            # write atomindices for each selection to list
            indices1 = []
            for at in sel1.atoms:
                indices1.append(at.index)
            indices2 = []
            for at in sel2.atoms:
                indices2.append(at.index)

            currentFrameContacts = []
            frame = ts.frame
            self.currentFrameNumber = ts.frame
            # print(frame)
            result = np.ndarray(shape=(len(sel1.positions), len(sel2.positions)), dtype=float)
            # distarray is the distance matrix between all atoms in sel1 and sel2
            # row = sel1, column = sel2
            distarray = distances.distance_array(sel1.positions, sel2.positions, box=None, result=result)
            contacts = np.where(distarray <= cutoff)
            # idx1 and idx2 correspond to a row,column in contacts, respectively
            # they do NOT correspond to a global atom index!
            for idx1, idx2 in itertools.izip(contacts[0], contacts[1]):
                convindex1 = indices1[idx1]  # idx1 converted to global atom indexing
                convindex2 = indices2[idx2]  # idx2 converted to global atom indexing
                # jump out of loop if hydrogen contacts are found, only contacts between heavy atoms are considered,
                # hydrogen bonds can still be detected!
                if re.match("H(.*)", self.name_array[convindex1]) or re.match("H(.*)", self.name_array[convindex2]):
                    continue
                    # distance between atom1 and atom2
                # check if residues are more than 4 apart, and in the same segment
                if selfInteraction:
                    if (self.resid_array[convindex1] - self.resid_array[convindex2]) < 5 and self.segids[convindex1] == self.segids[convindex2]:
                        continue
                distance = distarray[idx1, idx2]
                weight = self.weight_function(distance)

                # HydrogenBondAlgorithm
                hydrogenBonds = []
                # FF independent hydrogen bonds
                if (self.name_array[convindex1][0] in HydrogenBondAtoms.atoms and self.name_array[convindex2][0] in HydrogenBondAtoms.atoms):
                        # print("hbond? %s - %s" % (type_array[convindex1], type_array[convindex2]))
                        # search for hatom, check numbering in bond!!!!!!!!!!
                        b1 = self.bonds[convindex1]
                        b2 = self.bonds[convindex2]

                        # b1 = all_sel[convindex1].bonds
                        # b2 = all_sel[convindex2].bonds
                        # search for hydrogen atoms bound to atom 1
                        bondcount1 = 0
                        hydrogenAtomsBoundToAtom1 = []

                        # new code
                        for bnd in b1:
                            b = bnd.type
                            hydrogen = next((x for x in b if x.startswith("H")), 0)
                            # print(b)
                            if hydrogen != 0:
                                # print("h bond to atom1")
                                bondindices1 = b1.to_indices()[bondcount1]
                                # print bondindices1
                                # for j in bondindices1:
                                #     print(self.type_array[j+1])
                                hydrogenidx = next(
                                    (j for j in bondindices1 if self.name_array[j].startswith("H")), -1)
                                if hydrogenidx != -1:
                                    # print(self.type_array[hydrogenidx])
                                    hydrogenAtomsBoundToAtom1.append(hydrogenidx)
                            bondcount1 += 1
                        # search for hydrogen atoms bound to atom 2
                        bondcount2 = 0
                        hydrogenAtomsBoundToAtom2 = []
                        # print(b2)
                        for bnd2 in b2:
                            b = bnd2.type
                            hydrogen = next((x for x in b if x.startswith("H")), 0)
                            # print(b)
                            if hydrogen != 0:
                                # print("h bond to atom2")
                                bondindices2 = b2.to_indices()[bondcount2]
                                hydrogenidx = next(
                                    (k for k in bondindices2 if self.name_array[k].startswith("H")), -1)
                                if hydrogenidx != -1:
                                    # print(type_array[hydrogenidx])
                                    hydrogenAtomsBoundToAtom2.append(hydrogenidx)
                            bondcount2 += 1
                        # check hbond criteria for hydrogen atoms bound to first atom
                        for global_hatom in hydrogenAtomsBoundToAtom1:
                            conv_hatom = indices1.index(global_hatom)
                            # print(typeHeavy)
                            #
                            # TODO: FF independent version
                            # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[conv_hatom, idx2] <= hbondcutoff):
                            dist = distarray[conv_hatom, idx2]
                            if (dist <= hbondcutoff):
                                donorPosition = sel1.positions[idx1]
                                hydrogenPosition = sel1.positions[conv_hatom]
                                acceptorPosition = sel2.positions[idx2]
                                v1 = hydrogenPosition - acceptorPosition
                                v2 = hydrogenPosition - donorPosition
                                v1norm = np.linalg.norm(v1)
                                v2norm = np.linalg.norm(v2)
                                dot = np.dot(v1, v2)
                                angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                                # print(angle)
                                if angle >= hbondcutangle:
                                    # print("new hbond")
                                    new_hbond = HydrogenBond(convindex1, convindex2, global_hatom, dist, angle,
                                                             hbondcutoff,
                                                             hbondcutangle)
                                    hydrogenBonds.append(new_hbond)
                                    # print(str(convindex1) + " " + str(convindex2)
                                    # print("hbond found: %d,%d,%d"%(convindex1,global_hatom,convindex2))
                                    # print(angle)
                        for global_hatom in hydrogenAtomsBoundToAtom2:
                            conv_hatom = indices2.index(global_hatom)
                            # TODO: FF independent version
                            # if (typeHeavy == AtomHBondType.acc or typeHeavy == AtomHBondType.both) and (distarray[idx1, conv_hatom] <= hbondcutoff):
                            # FIXME: WTF?
                            # if (distarray[conv_hatom, idx2] <= hbondcutoff):
                            dist = distarray[idx1, conv_hatom]
                            if (dist <= hbondcutoff):
                                donorPosition = sel2.positions[idx2]
                                hydrogenPosition = sel2.positions[conv_hatom]
                                acceptorPosition = sel1.positions[idx1]
                                v1 = hydrogenPosition - acceptorPosition
                                v2 = hydrogenPosition - donorPosition
                                v1norm = np.linalg.norm(v1)
                                v2norm = np.linalg.norm(v2)
                                dot = np.dot(v1, v2)
                                angle = np.degrees(np.arccos(dot / (v1norm * v2norm)))
                                if angle >= hbondcutangle:
                                    new_hbond = HydrogenBond(convindex2, convindex1, global_hatom, dist, angle,
                                                             hbondcutoff,
                                                             hbondcutangle)
                                    hydrogenBonds.append(new_hbond)
                                    # print str(convindex1) + " " + str(convindex2)
                                    # print "hbond found: %d,%d,%d"%(convindex2,global_hatom,convindex1)
                                    # print angle
                                    # finalize
                newAtomContact = AtomContact(int(frame), float(distance), float(weight), int(convindex1),
                                             int(convindex2),
                                             hydrogenBonds)
                currentFrameContacts.append(newAtomContact)
            contactResults.append(currentFrameContacts)
        stop = time.time()

        # show trajectory information and selection information
        #print("trajectory with %d frames loaded" % len(u.trajectory))
        #print("Selection 1: ", len(sel1.positions), ", Selection2: ", len(sel2.positions))

        #print("analyzeTime: ", stop - start)
        # pickle.dump(contactResults, open("single_results_working.dat", "w"))
        # for f in contactResults:
            # print("working", len(f))
        print("distmatrix: ", stop-start)
        return contactResults

    def analyze_trackMolecule(self, contactResults, selindex, map):
        total = len(contactResults)
        counter = 1

        # list with list of tuples with contacts sorted by frame
        allSortedFrameContacts = []
        for frame in contactResults:
            currentFrameAcc = {}
            for cont in frame:
                if selindex == 1:
                    key1, key2 = self.makeKeyArraysFromMaps([], map, cont)
                    tit = Analyzer.make_single_title(key2)
                elif selindex == 2:
                    key1, key2 = self.makeKeyArraysFromMaps(map, [], cont)
                    tit = Analyzer.make_single_title(key1)
                if tit in currentFrameAcc:
                    currentFrameAcc[tit] += cont.weight
                else:
                    currentFrameAcc[tit] = 0
            sorted_frame_contacts = sorted(currentFrameAcc.items(), key=operator.itemgetter(1), reverse=1)
            allSortedFrameContacts.append(sorted_frame_contacts)
        return allSortedFrameContacts



    def analyze_contactResultsWithMaps(self, contactResults, map1, map2):
        """Analyzes contactsResults with the given maps."""

        #################################################
        # contactResults evaluation
        # only depending on map1, map2
        # part can be run without running the contact analysis algorithm again,
        # as it just prepares the results for displaying
        #################################################

        # Data structure to process:
        # contactResults (list of frames)
        # ---> frame (list of AtomContacts)
        # --------> AtomContact

        frame_contacts_accumulated = []
        # frame_contacts_accumulated (list of frames)
        # ---> frame_dict (dict)
        # --------> key vs. TempContactAccumulate

        # list of all contacts keys (= unique identifiers, determined by the given maps)
        start = time.time()
        allkeys = []
        total = len(contactResults)
        counter = 1
        for frame in contactResults:
            currentFrameAcc = {}
            for cont in frame:
                key1, key2 = self.makeKeyArraysFromMaps(map1, map2, cont)
                key = self.makeKeyFromKeyArrays(key1, key2)
                if key in currentFrameAcc:
                    currentFrameAcc[key].fscore += cont.weight
                    currentFrameAcc[key].contributingAtomContacts.append(cont)
                    if cont.idx1 in self.backbone:
                        currentFrameAcc[key].bb1score += cont.weight
                    else:
                        currentFrameAcc[key].sc1score += cont.weight
                    if cont.idx2 in self.backbone:
                        currentFrameAcc[key].bb2score += cont.weight
                    else:
                        currentFrameAcc[key].sc2score += cont.weight
                else:
                    currentFrameAcc[key] = TempContactAccumulate(key1, key2)
                    currentFrameAcc[key].fscore += cont.weight
                    currentFrameAcc[key].contributingAtomContacts.append(cont)
                    if cont.idx1 in self.backbone:
                        currentFrameAcc[key].bb1score += cont.weight
                    else:
                        currentFrameAcc[key].sc1score += cont.weight
                    if cont.idx2 in self.backbone:
                        currentFrameAcc[key].bb2score += cont.weight
                    else:
                        currentFrameAcc[key].sc2score += cont.weight
                if key not in allkeys:
                    allkeys.append(key)
            frame_contacts_accumulated.append(currentFrameAcc)
            self.frameUpdate.emit(float(counter) / float(total))
            counter += 1
        accumulatedContactsDict = {}
        stop = time.time()
        # print(stop - start)
        # accumulatedContactsDict (dict)
        # ---> key vs. list of TempContactAccumulated
        #
        # loop fills gaps with zero-score TempContactAccumulate of key if key is not occuring in a frame
        # provides clean data!
        start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if key not in frame_dict:  # puts empty score TempContactAccumulate in dict
                    key1, key2 = self.makeKeyArraysFromKey(key)
                    emptyCont = TempContactAccumulate(key1, key2)
                    emptyCont.fscore = 0
                    frame_dict[key] = emptyCont
                accumulatedContactsDict[key].append(frame_dict[key])

                # make a list of AccumulatedContacts from accumulatedContactsDict
        # probably, there is a much easier way to do that, but I am too tired at the moment and it works, though... (M)
        finalAccumulatedContacts = []  # list of AccumulatedContacts
        for key in accumulatedContactsDict:
            key1, key2 = self.makeKeyArraysFromKey(key)
            acc = AccumulatedContact(key1, key2)
            for tempContact in accumulatedContactsDict[key]:
                acc.addScore(tempContact.fscore)
                acc.addContributingAtoms(tempContact.contributingAtomContacts)
                acc.bb1 += tempContact.bb1score
                acc.bb2 += tempContact.bb2score
                acc.sc1 += tempContact.sc1score
                acc.sc2 += tempContact.sc2score
            finalAccumulatedContacts.append(acc)
            # print(key, acc.bb1, acc.bb2, acc.sc1, acc.sc2)
            # print(len(acc.scoreArray))
        stop = time.time()
        # print(stop - start)
        return finalAccumulatedContacts

    def analysisEventListener(self):
        """Event listener for the progress bar in MainWindow."""
        while self.analysis_state:
            progress = 0
            for each in analysisProgressDict.keys():
                progress += analysisProgressDict[each]
            progress = float(progress) / float(self.totalFramesToProcess)
            if progress > 0:
                self.frameUpdate.emit(progress)

            if progress == 1.0:
                for each in analysisProgressDict.keys():
                    analysisProgressDict[each] = 0
                progress = 0
                self.analysis_state = False

    def analyze_contactResultsWithMaps_Parallel(self, map1, map2, nproc):
        """Analyzes contactsResults with the given maps, using nproc threads."""
        start = time.time()
        self.totalFramesToProcess = len(self.contactResults)
        results = []
        rank = 0
        manager = multiprocessing.Manager()
        d = manager.list(self.getTrajectoryData())
        all_chunk = chunks(self.contactResults, nproc)
        pool = LoggingPool(nproc)
        print("Running on %d cores" % nproc)
        for c in all_chunk:
            results.append(pool.apply_async(loop_frame, args=(c, map1, map2, d, rank)))
            rank += 1
        # self.totalFramesToProcess = len(contResults)
        self.analysis_state = True
        self.analysisEventListener()
        pool.close()
        pool.join()
        self.analysis_state = False
        stop = time.time()
        # print("time: ", str(stop - start), rank)
        # print(str(len(c)), rank)
        allkeys = []
        frame_contacts_accumulated = []
        # print(len(results))
        for res in results:
            rn = res.get()
            allkeys.extend(rn[0])
            frame_contacts_accumulated.extend(rn[1])
        accumulatedContactsDict = {}
        start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if key not in frame_dict:  # puts empty score TempContactAccumulate in dict
                    key1, key2 = makeKeyArraysFromKey(key)
                    emptyCont = TempContactAccumulate(key1, key2)
                    emptyCont.fscore = 0
                    frame_dict[key] = emptyCont
                accumulatedContactsDict[key].append(frame_dict[key])
        finalAccumulatedContacts = []  # list of AccumulatedContacts
        for key in accumulatedContactsDict:
            key1, key2 = makeKeyArraysFromKey(key)
            acc = AccumulatedContact(key1, key2)
            for tempContact in accumulatedContactsDict[key]:
                acc.addScore(tempContact.fscore)
                acc.addContributingAtoms(tempContact.contributingAtomContacts)
                acc.bb1 += tempContact.bb1score
                acc.bb2 += tempContact.bb2score
                acc.sc1 += tempContact.sc1score
                acc.sc2 += tempContact.sc2score
            finalAccumulatedContacts.append(acc)
        # stop = time.time()
        # print(stop - start)
        glob_stop = time.time()
        # print(glob_stop - start)
        # self.progressWidget.hide()
        return finalAccumulatedContacts
