from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np

from .Biochemistry import *


class Operator(object):
    """Defines a comparison operator for contact filtering."""
    greater, smaller, equal, nequal = range(4)
    mapping = {u"greater": greater, u"smaller": smaller, u"equal": equal, u"not equal": nequal}

    def compare(self, value1, value2, operator):
        if operator == self.greater:
            return value1 > value2
        elif operator == self.smaller:
            return value1 < value2
        elif operator == self.equal:
            return value1 == value2
        elif operator == self.nequal:
            return value1 != value2


class FrameFilter(object):
    """Filters contact with a given frame range."""
    def __init__(self, name):
        self.name = name

    @staticmethod
    def extractFrameRange(contacts, frameRange):
        lower = frameRange[0]
        upper = frameRange[1]
        for c in contacts:
            newScores = c.scoreArray[lower:upper]
            newAtoms = c.contributingAtoms[lower:upper]
            c.contributingAtoms = newAtoms
            c.scoreArray = newScores
        return contacts


class NameFilter(object):
    """Filters contacts by name."""

    def __init__(self, name):
        self.name = name

    @staticmethod
    def filterContactsByName(contacts, nameA, nameB, mapindex):
        """Returns contacts, when their key is included in nameA and nameB."""
        filtered = []
        for c in contacts:
            add = False
            # TODO: probably replace with if/else
            try:
                prop1 = c.key1[mapindex]
                prop2 = c.key2[mapindex]
            except IndexError:
                filtered.append(c)
                continue
            if nameA.lower() != u'all' and nameB.lower() != u'all':
                splitA = nameA.split(u",")
                splitB = nameB.split(u",")
                if prop1 in splitA and prop2 in splitB:
                    add = True
            elif nameA.lower() != u'all' and nameB.lower() == u'all':
                splitA = nameA.split(u",")
                if prop1 in splitA:
                    add = True
            elif nameA.lower() == u'all' and nameB.lower() != u'all':
                splitB = nameB.split(u",")
                if prop2 in splitB:
                    add = True
            else:
                add = True

            if add:
                filtered.append(c)
        return filtered


class RangeFilter(object):
    """Filters contacts with a given range."""
    def __init__(self, name):
        self.name = name

    @staticmethod
    def numberInRanges(number, ranges):
        result = False
        for r in ranges:
            if number in r:
                result = True
        return result

    def filterByRange(self, contacts, residRangeA, residRangeB, mapindex):
        """Returns contacts that lie in residRangeA and residRangeB."""
        splitA = residRangeA.split(u",")
        splitB = residRangeB.split(u",")

        notAllA = residRangeA.lower() != u'all'
        notAllB = residRangeB.lower() != u'all'

        if notAllA:
            aRanges = []
            for ran in splitA:
                r = ran.split(u"-")
                if len(r) == 1:
                    aRanges.append(range(int(r[0]), int(r[0])+1))
                else:
                    aRanges.append(range(int(r[0]), int(r[1]) + 1))
        if notAllB:
            bRanges = []
            for ran in splitB:
                r = ran.split(u"-")
                print(r)
                if len(r) == 1:
                    bRanges.append(range(int(r[0]), int(r[0])+1))
                else:
                    bRanges.append(range(int(r[0]), int(r[1]) + 1))

        filtered = []
        for c in contacts:
            # TODO: probably replace with if/else
            try:
                prop1 = int(c.key1[mapindex])
                prop2 = int(c.key2[mapindex])
            except ValueError:
                filtered.append(c)
                continue
            add = False
            if notAllA and notAllB:
                if self.numberInRanges(prop1, aRanges) and self.numberInRanges(prop2, bRanges):
                    add = True
            elif notAllA and not notAllB:
                if self.numberInRanges(prop1, aRanges):
                    add = True
            elif not notAllA and notAllB:
                if self.numberInRanges(prop2, bRanges):
                    add = True
            else:
                add = True

            if add:
                filtered.append(c)
        return filtered


class BinaryFilter(object):
    """Implements a binary filter with a given operator."""
    def __init__(self, name, operator, value):
        self.name = name
        self.operator = Operator.mapping[operator]
        self.value = value

    def filterContacts(self, contacts):
        pass


class OnlyFilter(object):
    """Implements a filter, that only selects contacts with specific properties, e.g.: hydrophobic."""
    def __init__(self, name, operator, value):
        self.name = name
        self.operator = operator
        self.value = value

    def filterContacts(self, contacts):
        filtered = []
        if self.operator == "hbonds":
            for c in contacts:
                hb = c.hbondFramesScan()
                for bla in hb:
                    if bla > 0:
                        filtered.append(c)
                        break
        elif self.operator == "hydrophobic":
            for c in contacts:
                if c.determine_ctype() == ContactType.hydrophobic:
                    filtered.append(c)
        elif self.operator == "saltbridges":
            for c in contacts:
                if c.determine_ctype() == ContactType.saltbr:
                    filtered.append(c)
        elif self.operator == "other":
            for c in contacts:
                if c.determine_ctype() == ContactType.other:
                    filtered.append(c)
        return filtered


class TotalTimeFilter(BinaryFilter):
    """Binary filter concerning the total contact time."""
    def __init__(self, name, operator, value):
        super(TotalTimeFilter, self).__init__(name, operator, value)

    def filterContacts(self, contacts):
        filtered = []
        op = Operator()
        for c in contacts:
            if op.compare(c.total_time(1, 0), self.value, self.operator):
                filtered.append(c)
        print(unicode(len(filtered)))
        return filtered


class ScoreFilter(BinaryFilter):
    """Compares contact score of every frame, only adds contact if true for all frames."""
    def __init__(self, name, operator, value, ftype):
        super(ScoreFilter, self).__init__(name, operator, value)
        self.ftype = ftype

    def filterContacts(self, contacts):
        filtered = []
        op = Operator()
        if self.ftype == u"Mean":
            for c in contacts:
                mean = c.mean_score()
                if op.compare(mean, self.value, self.operator):
                    filtered.append(c)
        elif self.ftype == u"Median":
            for c in contacts:
                med = c.median_score()
                if op.compare(med, self.value, self.operator):
                    filtered.append(c)
        elif self.ftype == u"HB %":
            for c in contacts:
                med = c.hbond_percentage()
                print(med)
                if op.compare(med, self.value, self.operator):
                    filtered.append(c)
        return filtered


class SortingOrder(object):
    """Defines the sorting order, ascending or descending."""
    ascending, descending = range(2)
    mapping = {u"asc.": ascending, u"desc.": descending}


class Sorting(object):
    """Performs the sorting with respect to a specific key, e.g. mean lifetime."""
    def __init__(self, name, key, descending):
        self.name = name
        self. key = key
        self.descending = descending
        self.threshold = 0
        self.nspf = 0

    def setThresholdAndNsPerFrame(self, threshold, nspf):
        self.threshold = threshold
        self.nspf = nspf

    def sortContacts(self, contacts):
        sortedContacts = []
        if self.key == u"mean":
            sortedContacts = sorted(contacts, key=lambda c: c.meanScore, reverse=self.descending)
        elif self.key == u"median":
            sortedContacts = sorted(contacts, key=lambda c: c.medianScore, reverse=self.descending)
        elif self.key == u"bb/sc type":
            sortedContacts = sorted(contacts, key=lambda c: c.backboneSideChainType, reverse=self.descending)
        elif self.key == u"contact type":
            sortedContacts = sorted(contacts, key=lambda c: c.contactType, reverse=self.descending)
        elif self.key == u"resid 1":
            try:
                prop1 = contacts[0].key1[AccumulationMapIndex.resid]
                if prop1 != "none":
                    sortedContacts = sorted(contacts, key=lambda c: int(c.key1[AccumulationMapIndex.resid]), reverse=self.descending)
                else:
                    sortedContacts = contacts
            except IndexError:
                pass
        elif self.key == u"resid 2":
            try:
                prop2 = contacts[0].key2[AccumulationMapIndex.resid]
                if prop2 != "none":
                    sortedContacts = sorted(contacts, key=lambda c: int(c.key2[AccumulationMapIndex.resid]), reverse=self.descending)
                else:
                    sortedContacts = contacts
            except IndexError:
                pass
        elif self.key == u"total time":
            for con in contacts:
                con.total_time(self.nspf, self.threshold)
            sortedContacts = sorted(contacts, key=lambda c: c.ttime, reverse=self.descending)
        elif self.key == u"mean lifetime":
            for con in contacts:
                con.mean_life_time(self.nspf, self.threshold)
            sortedContacts = sorted(contacts, key=lambda c: c.meanLifeTime, reverse=self.descending)
        elif self.key == u"median lifetime":
            for con in contacts:
                con.median_life_time(self.nspf, self.threshold)
            sortedContacts = sorted(contacts, key=lambda c: c.medianLifeTime, reverse=self.descending)
        return sortedContacts
