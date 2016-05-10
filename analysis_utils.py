#### utilities
### should be moved to a separate file in the future project...

# utility for memory measurement
import numpy as np

suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']


def humansize(nbytes):
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes) - 1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])


# find a string in s between the strings first and last
def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""

# weight function to score contact distances
def weight_function(value):
    return (1.0) / (1.0 + np.exp(5.0 * (value - 4.0)))


## enum and mapping for atom properties
# used to dynamically define keys and have a bijective nomenclature for all properties
class AccumulationMapIndex():
    index, atype, name, resid, resname, segid = range(6)
    mapping = ["i.", "t.", "nm.", "r.", "rn.", "s."]


## maps contain information wether to consider an atom's field for contact accumulation
# map1 and map2 contain six boolean values each, cf. AccumulationMapIndex
# for a given contact, the corresponding value to a field is written to keys1 and keys2, respectively
# example input:
# map1 = [0,0,0,1,1,0]
# map2 = [0,0,0,1,1,0], meaning that residue and resname should be used for contact accumulation
# contact: idx1,idx2
# results: (example!)
#	keys1=["none","none","none","14", "VAL", "none"]
#	keys2=["none","none","none","22", "ILE, "none"]
def makeKeyArraysFromMaps(map1, map2, contact):
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


# convert a key back to two key arrays
# cf. comments on makeKeyFromKeyArrays and makeKeyArraysFromMaps
# "inverse" function of makeKeyFromKeyArrays
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


## input two key arrays as explained above
#	example:
#	keys1=["none","none","none","14", "VAL", "none"]
#	keys2=["none","none","none","22", "ILE, "none"]
#	returns a human readable key with the mapping identifiers in AccumulationMapIndex
#	in the given example data:
#	key="r.14rn.VAL-r.22rn.ILE"
#	key is used to accumulated contacts in a dictionary (= a contact's unique identifier)
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