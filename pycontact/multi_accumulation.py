import numpy as np
import MDAnalysis
from mdanalysis import *
import math
import time
import multiprocessing

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
    print 'RATIO', ratio
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

def loop_frame(contacts,map1,map2,trajArgs):
    allkeys = []
    results = []
    global backbone
    global type_array
    global name_array
    global resid_array
    global resname_array
    global segids
    backbone,type_array,name_array,resid_array,resname_array,segids = trajArgs[-1],trajArgs[3],trajArgs[2],trajArgs[1],trajArgs[0],trajArgs[4]
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
    return [allkeys,results]


# if __name__ == '__main__':
#     glob_start = time.time()
#     start = time.time()
#     import pickle
#     # importDict = pickle.load(open("defaultsession", "rb"))
#     importDict = pickle.load(open("/Users/maximilianscheurer/Dropbox/TCBG/ba/data/prot_ubp6/yeast_proteasome_ubp6_session", "rb"))
#     contResults = importDict["analyzer"][-1]
#     trajArgs = importDict["trajectory"]
#     nproc = 8
#     tasks = []
#     results = []
#     rank = 0
#     map1=[1,1,1,1,1,1]
#     map2=[1,1,1,1,1,1]
#     # backbone,type_array,name_array,resid_array,resname_array,segids = trajArgs[-1],trajArgs[3],trajArgs[2],trajArgs[1],trajArgs[0],trajArgs[4]
#     manager = multiprocessing.Manager()
#     d=manager.list(trajArgs)
#     all_chunk = chunks(contResults,nproc)
#     pool = multiprocessing.Pool(nproc)
#     for c in all_chunk:
#         results.append( pool.apply_async( loop_frame, args=(c,map1,map2,d)) )
#         rank +=1
#     pool.close()
#     pool.join()
#     stop = time.time()
#     print "time: ", str(stop-start), rank
#     print str(len(c)), rank
#     allkeys = []
#     frame_contacts_accumulated = []
#     print len(results)
#     for res in results:
#         rn = res.get()
#         allkeys.extend(rn[0])
#         frame_contacts_accumulated.extend(rn[1])
#     accumulatedContactsDict = {}
#     #   start = time.time()
#     for key in allkeys:
#         accumulatedContactsDict[key] = []
#         for frame_dict in frame_contacts_accumulated:
#             if not key in frame_dict:  # puts empty score TempContactAccumulate in dict
#                 key1, key2 = makeKeyArraysFromKey(key)
#                 emptyCont = TempContactAccumulate(key1, key2)
#                 emptyCont.fscore = 0
#                 frame_dict[key] = emptyCont
#             accumulatedContactsDict[key].append(frame_dict[key])
#     finalAccumulatedContacts = []  # list of AccumulatedContacts
#     for key in accumulatedContactsDict:
#         key1, key2 = makeKeyArraysFromKey(key)
#         acc = AccumulatedContact(key1, key2)
#         for tempContact in accumulatedContactsDict[key]:
#             acc.addScore(tempContact.fscore)
#             acc.addContributingAtoms(tempContact.contributingAtomContacts)
#             acc.bb1 += tempContact.bb1score
#             acc.bb2 += tempContact.bb2score
#             acc.sc1 += tempContact.sc1score
#             acc.sc2 += tempContact.sc2score
#         finalAccumulatedContacts.append(acc)
#     # stop = time.time()
#     # print stop - start
#     glob_stop = time.time()
#     print glob_stop - glob_start
