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
