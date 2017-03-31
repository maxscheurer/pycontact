"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
matplotlib.use("qt5agg")

from PyContact.core.Scripting import PyContactJob, JobConfig
from PyContact.core.Biochemistry import AccumulationMapIndex

jobAllMD = PyContactJob("./PyContact/exampleData/rpn11_ubq_interface-ionized.psf",
                   "./PyContact/exampleData/short.dcd", "md_all", JobConfig(5.0, 2.5, 120, [0,0,0,1,1,0], [0,0,0,1,1,0], "segid UBQ", "segid RN11"))

jobAllMD.runJob(1)
contacts = jobAllMD.analyzer.finalAccumulatedContacts

map1 = [0,0,0,1,1,0]
map2 = [0,0,0,1,1,0]
label1 = "UBQ"
label2 = "RN11"
attribute = "test"
threshold = 0.0
nsPerFrame = 1.0

minmaxresids1 = []
minmaxresids2 = []

for cont in contacts:
    minmaxresids1.append(int(cont.key1[AccumulationMapIndex.resid]))
    minmaxresids2.append(int(cont.key2[AccumulationMapIndex.resid]))

x = np.arange(np.min(minmaxresids1), np.max(minmaxresids1) + 1)
y = np.arange(np.min(minmaxresids2), np.max(minmaxresids2) + 1)

fig, ax = plt.subplots()

data = np.ones((len(x), len(y)))
# data = np.random.rand(len(x), len(y))
print(data.shape)
print(data)

# cax = ax.matshow(data, cmap=plt.cm.gray)
cax = None
ttl = ax.set_title("test")
plt.draw()
minx = np.min(minmaxresids1)
miny = np.min(minmaxresids2)
rng = np.arange(len(contacts[0].scoreArray))
# ttl = ax.text(30, -10, ' ')


def animate(i):
    print(i)
    global cax
    data = np.zeros((len(x), len(y)))
    attribute = "Mean Score"
    if attribute == "Mean Score":
        for c in contacts:
            r1 = int(c.key1[AccumulationMapIndex.resid]) - minx
            r2 = int(c.key2[AccumulationMapIndex.resid]) - miny
            data[r1, r2] = c.scoreArray[i]
            print(data[r1, r2])
    cax.set_data(data)  # update the data
    ttl.set_text(str(i))
    plt.draw()
    return cax,


# Init only required for blitting to give a clean slate.
def init():
    global cax
    data = np.zeros((len(x), len(y)))

    frame = 0
    print("frame: ", frame)
    attribute = "Mean Score"
    if attribute == "Mean Score":
        for c in contacts:
            r1 = int(c.key1[AccumulationMapIndex.resid]) - minx
            r2 = int(c.key2[AccumulationMapIndex.resid]) - miny
            data[r1, r2] = c.scoreArray[frame]
            print(data[r1, r2])

    # data.fill(10)
    # data = np.zeros((len(x), len(y)))
    cax = ax.matshow(data, cmap=plt.cm.gray)
    ttl.set_text("init")
    plt.draw()
    return cax,

ani = animation.FuncAnimation(fig, animate, rng, init_func=init, blit=True, repeat=False)
plt.show()
