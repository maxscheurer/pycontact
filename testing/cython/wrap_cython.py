# testing wrapper

import matplotlib.pyplot as plt
from wrap_cy import bla, test_sasa

a = bla(2)
print(a)

import time
start = time.time()
r = test_sasa()
stop = time.time()
print(stop-start)
plt.plot(r)
plt.show()
