import numpy as np
from biochemistry import *

class Operator:
    greater, smaller, equal, nequal = range(4)
    mapping = {"greater": greater, "smaller": smaller, "equal": equal, "not equal": nequal}

    def compare(self, value1, value2, operator):
        if operator == self.greater:
            return (value1 > value2)
        elif operator == self.smaller:
            return (value1 < value2)
        elif operator == self.equal:
            return (value1 == value2)
        elif operator == self.nequal:
            return (value1 != value2)


class BinaryFilter:
    def __init__(self, name, operator, value):
        self.name = name
        self.operator = Operator.mapping[operator]
        self.value = value

    def filterContacts(self,contacts):
        pass


class TotalTimeFilter(BinaryFilter):
    def __init__(self,name, operator, value):
        super(TotalTimeFilter, self).__init__(name, operator, value)

    def filterContacts(self,contacts):
        filtered = []
        op = Operator()
        for c in contacts:
            if op.compare(c.total_time(1, 0), self.value, self.operator):
                filtered.append(c)
        print(str(len(filtered)))
        return filtered

# filter compares contact score of every frame, only adds contact if true for all frames
class ScoreFilter(BinaryFilter):
    def __init__(self, name, operator, value, ftype):
        super(ScoreFilter, self).__init__(name, operator, value)
        self.ftype = ftype

    def filterContacts(self, contacts):
        filtered = []
        op = Operator()
        if self.ftype == "mean":
            for c in contacts:
                mean = c.mean_score()
                if op.compare(mean, self.value, self.operator):
                    filtered.append(c)
        elif self.ftype == "median":
            for c in contacts:
                med = c.median_score()
                if op.compare(med, self.value, self.operator):
                    filtered.append(c)
        return filtered

class SortingOrder:
    ascending, descending = range(2)
    mapping = {"asc": ascending, "desc": descending}

class Sorting:
    def __init__(self, name, key, descending):
        self.name = name
        self. key = key
        self.descending = descending

    def setThresholdAndNsPerFrame(self,threshold,nspf):
        self.threshold = threshold
        self.nspf = nspf

    def sortContacts(self, contacts):
        sortedContacts = []
        if self.key == "mean":
            sortedContacts = sorted(contacts, key=lambda c: c.meanScore, reverse=self.descending)
        elif self.key == "median":
            sortedContacts = sorted(contacts, key=lambda c: c.medianScore, reverse=self.descending)
        elif self.key == "bb/sc type":
            sortedContacts = sorted(contacts, key=lambda c: c.backboneSideChainType, reverse=self.descending)
        elif self.key == "total time":
            for con in contacts:
                con.total_time(self.nspf,self.threshold)
            sortedContacts = sorted(contacts, key=lambda c: c.ttime, reverse=self.descending)
        return sortedContacts