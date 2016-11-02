#usage: mpirun -np 4 python mpi_first.py
from mpi4py import MPI
import numpy as np
import MDAnalysis
from mdanalysis import *
import math
import time
glob_start = time.time()

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
	global type_array,name_array,resid_array,resname_array,segids
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

def loop_frame(contacts,map1,map2):
	global backbone
	allkeys = []
	results = []
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

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
	start = time.time()
	import pickle
	importDict = pickle.load(open("defaultsession", "rb"))
	# importDict = pickle.load(open("/Users/maximilianscheurer/Dropbox/TCBG/ba/data/prot_ubp6/yeast_proteasome_ubp6_session", "rb"))
	stop = time.time()
	print stop-start
	contResults = importDict["analyzer"][-1]
	trajArgs = importDict["trajectory"]
	all_chunk = list(chunks(contResults,size))
	for rk in range(1,size):
		comm.send(trajArgs, dest=rk,tag=(math.pow(rk,2)+7))
		# comm.send(all_chunk[rk], dest=rk,tag=(math.pow(rk,2)+9))
else:
	all_chunk = None
	rec = time.time()
	trajArgs = comm.recv(source=0,tag=(math.pow(rank,2)+7))
	rec_end = time.time()
	print 'rec',(rec_end-rec),rank
	# all_chunk = comm.recv(source=0,tag=(math.pow(rank,2)+9))
	# [self.resname_array,self.resid_array,self.name_array,self.type_array,self.segids,self.backbone]
	# map1=[1,1,1,1,1]
	# map2=[1,1,1,1,1]
	# start = time.time()
	# results = loop_frame(all_chunk,map1,map2,trajArgs[-1],trajArgs[3],trajArgs[2],trajArgs[1],trajArgs[0],trajArgs[4])
	# stop = time.time()
	# print "time: ", str(stop-start), rank
	# print str(len(all_chunk)), rank
all_chunk = comm.scatter(all_chunk,root=0)
map1=[0,0,0,1,1,0]
map2=[0,0,0,1,1,0]
start = time.time()
backbone,type_array,name_array,resid_array,resname_array,segids = trajArgs[-1],trajArgs[3],trajArgs[2],trajArgs[1],trajArgs[0],trajArgs[4]
results = loop_frame(all_chunk,map1,map2)
stop = time.time()
print "time: ", str(stop-start), rank
print str(len(all_chunk)), rank

results = comm.gather(results,root=0)

if rank == 0:
	allkeys = []
	frame_contacts_accumulated = []
	print len(results)
	for rn in results:
		allkeys.extend(rn[0])
		frame_contacts_accumulated.extend(rn[1])
	accumulatedContactsDict = {}
	start = time.time()
        for key in allkeys:
            accumulatedContactsDict[key] = []
            for frame_dict in frame_contacts_accumulated:
                if not key in frame_dict:  # puts empty score TempContactAccumulate in dict
                    key1, key2 = makeKeyArraysFromKey(key)
                    emptyCont = TempContactAccumulate(key1, key2)
                    emptyCont.fscore = 0
                    frame_dict[key] = emptyCont
                accumulatedContactsDict[key].append(frame_dict[key])

                # make a list of AccumulatedContacts from accumulatedContactsDict
        # probably, there is a much easier way to do that, but I am too tired at the moment and it works, though... (M)
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
            # print key, acc.bb1, acc.bb2, acc.sc1, acc.sc2
            # print len(acc.scoreArray)
        stop = time.time()
        print stop - start
        glob_stop = time.time()
        print glob_stop - glob_start
# all_chunk = comm.scatter(all_chunk, root=0)
#print 'rank',rank,'has data:',data
# newData = comm.gather(arguments,root=0)

# if rank == 0:
   # print 'master:',newData
