from. Biochemistry import AccumulationMapIndex


def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""

class KeyManager:
    """docstring for [object Object]."""
    def __init__(self, name_array, resid_array, resname_array, segids):
        self.name_array = name_array
        self.resid_array = resid_array
        self.resname_array = resname_array
        self.segids = segids

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
    def makeKeyArraysFromKey(key):
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
