import pickle

s = pickle.load(open("single_results.dat"))
print len(s)
count = 0
for element in s:
    for c in element:
        count += len(c.hbondinfo)
        if len(c.hbondinfo) > 0:
            c.hbondinfo[0].toString()
print count


p = pickle.load(open("parallel_results.dat"))
print len(s)
count = 0
for element in p:
    for c in element:
        count += len(c.hbondinfo)
        if len(c.hbondinfo) > 0:
            c.hbondinfo[0].toString()

print count
