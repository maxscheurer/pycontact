from __future__ import print_function
import math
import time
import multiprocessing

import numpy as np
import MDAnalysis
from MDAnalysis.analysis import distances

from .Biochemistry import *
analysisProgressManager = multiprocessing.Manager()
analysisProgressDict = analysisProgressManager.dict()

def find_between(s, first, last):
        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

def makeKeyArraysFromKey(key):
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
            currentValue = find_between(keystring1, current, nextCurrent)
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
            currentValue = find_between(keystring2, current, nextCurrent)
            if currentValue == "":
                key2.append("none")
            else:
                key2.append(currentValue)
    return [key1, key2]

def chunks_old(l, n):
    """Yield successive n-sized chunks from l."""
    ratio = int(math.floor(len(l)/n))
    print("RATIO", ratio)
    current_chunk = 1
    for i in range(0, len(l), ratio):
        if current_chunk == n:
            yield l[i:-1]
        else:
            yield l[i:i+ratio]
        current_chunk += 1

def chunks(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

def makeKeyArraysFromMaps(map1, map2, contact):
    # global type_array,name_array,resid_array,resname_array,segids
    idx1 = contact.idx1
    idx2 = contact.idx2
    counter = 0
    keys1 = []
    for val in map1:
        if val == 1:
            if counter == AccumulationMapIndex.index:
                keys1.append(idx1)
            elif counter == AccumulationMapIndex.atype:
                keys1.append(type_array[idx1])
            elif counter == AccumulationMapIndex.name:
                keys1.append(name_array[idx1])
            elif counter == AccumulationMapIndex.resid:
                keys1.append(resid_array[idx1])
            elif counter == AccumulationMapIndex.resname:
                keys1.append(resname_array[idx1])
            elif counter == AccumulationMapIndex.segid:
                keys1.append(segids[idx1])
        else:
            keys1.append("none")
        counter += 1
    counter = 0
    keys2 = []
    for val in map2:
        if val == 1:
            if counter == AccumulationMapIndex.index:
                keys2.append(idx2)
            elif counter == AccumulationMapIndex.atype:
                keys2.append(type_array[idx2])
            elif counter == AccumulationMapIndex.name:
                keys2.append(name_array[idx2])
            elif counter == AccumulationMapIndex.resid:
                keys2.append(resid_array[idx2])
            elif counter == AccumulationMapIndex.resname:
                keys2.append(resname_array[idx2])
            elif counter == AccumulationMapIndex.segid:
                keys2.append(segids[idx2])
        else:
            keys2.append("none")
        counter += 1
    return [keys1, keys2]

def makeKeyFromKeyArrays(key1, key2):
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

def loop_frame(contacts,map1,map2,trajArgs,rank):
    allkeys = []
    results = []
    global backbone
    global type_array
    global name_array
    global resid_array
    global resname_array
    global segids
    frames_processed = 0
    backbone,type_array,name_array,resid_array,resname_array,segids = trajArgs[5],trajArgs[3],trajArgs[2],trajArgs[1],trajArgs[0],trajArgs[4]
    for frame in contacts:
        currentFrameAcc = {}
        for cont in frame:
            key1, key2 = makeKeyArraysFromMaps(map1, map2, cont)
            key = makeKeyFromKeyArrays(key1, key2)
            # global rank
            # print key, rank
            if key in currentFrameAcc:
                currentFrameAcc[key].fscore += cont.weight
                currentFrameAcc[key].contributingAtomContacts.append(cont)
                if cont.idx1 in backbone:
                    currentFrameAcc[key].bb1score += cont.weight
                else:
                    currentFrameAcc[key].sc1score += cont.weight
                if cont.idx2 in backbone:
                    currentFrameAcc[key].bb2score += cont.weight
                else:
                    currentFrameAcc[key].sc2score += cont.weight
            else:
                currentFrameAcc[key] = TempContactAccumulate(key1, key2)
                currentFrameAcc[key].fscore += cont.weight
                currentFrameAcc[key].contributingAtomContacts.append(cont)
                if cont.idx1 in backbone:
                    currentFrameAcc[key].bb1score += cont.weight
                else:
                    currentFrameAcc[key].sc1score += cont.weight
                if cont.idx2 in backbone:
                    currentFrameAcc[key].bb2score += cont.weight
                else:
                    currentFrameAcc[key].sc2score += cont.weight
            if not key in allkeys:
                allkeys.append(key)
        results.append(currentFrameAcc)
        frames_processed += 1
        analysisProgressDict[rank] = frames_processed
    return [allkeys,results]
