class ContactType:
    hbond, saltbr, hydrophobic = range(3)

class ResidueType:
    positive, negative, unpolar, polar = range(4)



class Contact:
    def __init__(self, resA, residA, resB, residB, scoreArray):
        self.resA = resA
        self.resB = resB
        self.residA = residA
        self.residB = residB
        self.scoreArray = scoreArray
        self.title = self.resA + self.residA + "-" + self.resB + self.residB
        self.type = determine_ctype(self.resA, self.resB)

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
    return 0