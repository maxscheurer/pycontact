import collections
compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

class ResidueType:
    positive, negative, nonpolar, polar, other = range(5)
    mapping = {"asp": negative, "arg": positive, "lys": positive, "glu": negative}

class ContactType:
    saltbr, hydrophobic, other = range(3)
    mapping = [[ResidueType.positive, ResidueType.negative],[ResidueType.nonpolar, ResidueType.nonpolar]]
    colors = ["rgb(255, 0,0)", "rgb(0, 0,255)", "rgb(255, 255 ,255)"]

class Residue:
    def __init__(self, name):
        self.name = name.lower()
        if self.name in ResidueType.mapping:
            self.type = ResidueType.mapping[self.name]
        else:
            self.type = ResidueType.other

class Contact:
    def __init__(self, resA, residA, resB, residB, scoreArray):
        self.resA = resA
        self.resB = resB
        self.residA = residA
        self.residB = residB
        self.scoreArray = scoreArray
        self.title = self.resA + self.residA + "-" + self.resB + self.residB
        self.residueA = Residue(self.resA)
        self.residueB = Residue(self.resB)
        self.type = determine_ctype(self.residueA, self.residueB)

    def framenumber(self):
        return len(self.scoreArray)

    def total_time(self, ns_per_frame, threshold):
        time = 0
        for score in self.scoreArray:
            if score > threshold:
                time += ns_per_frame
        self.ttime = time
        return self.ttime

def determine_ctype(resA, resB):
    if compare([resA.type,resB.type],ContactType.mapping[ContactType.saltbr]):
        return ContactType.saltbr
    else:
        return ContactType.other

