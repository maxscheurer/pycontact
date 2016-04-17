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
        frame_info = np.array(single_elements[4:], dtype=float)
        newContact = Contact(resA, residA, resB, residB, frame_info)
        contacts.append(newContact)
    return contacts