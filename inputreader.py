import numpy as np
from biochemistry import *

def makeContactFromLines(lines):
    contacts = []
    for l in lines:
        single_elements = l.split()
        resA = single_elements[0]
        residA = single_elements[1]
        resB = single_elements[2]
        residB = single_elements[3]
        bb1 = single_elements[4]
        sc1 =single_elements[5]
        bb2 = single_elements[6]
        sc2 = single_elements[7]
        frame_info = np.array(single_elements[8:], dtype=float)
        newContact = Contact(resA, residA, resB, residB, bb1, sc1, bb2, sc2, frame_info)
        contacts.append(newContact)
    return contacts