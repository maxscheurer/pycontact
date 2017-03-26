from __future__ import print_function
from PyContact.core.Scripting import (PyContactJob, JobConfig)
from PyContact.core.ContactFilters import ScoreFilter
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

jobAllMD = PyContactJob("/home/max/Projects/pycontact/PyContact/exampleData/rpn11_ubq_interface-ionized.psf","/home/max/Projects/pycontact/PyContact/exampleData/short.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,0,1,1,0], [0,0,0,1,1,0], "segid UBQ", "segid RN11"))

jobAllMD.runJob(1)
# jobAllMD.writeSessionToFile()
contacts = jobAllMD.analyzer.finalAccumulatedContacts

filter = ScoreFilter("score", "greater", 0.0, u"Median")
contacts = filter.filterContacts(contacts)

N = len(contacts)
frames = len(jobAllMD.analyzer.contactResults)
print("frms:", frames)
covMatrix = np.zeros((N, N))


idx1 = 0
r_vec = np.zeros((N, frames))
r_names = np.zeros(N).tolist()
for c1 in contacts:
    idx2 = 0
    for c2 in contacts:
        for f in range(frames):
            covMatrix[idx1, idx2] += (c1.scoreArray[f] - c1.mean_score()) * (c2.scoreArray[f] - c2.mean_score())
            r_vec[idx1, f] = c1.scoreArray[f]
            r_names[idx1] = c1.human_readable_title()
        covMatrix[idx1, idx2] /= frames
        idx2 += 1
    idx1 += 1

# plt.matshow(covMatrix)
# plt.show()

pca_evalues, pca_evecs = LA.eigh(covMatrix, UPLO='L')
idx = pca_evalues.argsort()[::-1]
pca_evalues_sorted = pca_evalues[idx]
pca_evecs_sorted = pca_evecs[:, idx]

percentages = pca_evalues_sorted / (np.sum(pca_evalues_sorted))
# print(percentages)

pc_frames = []
pcOfInterest = 0

# ???? auto-correlation
for p in range(3):
    pc_list = []
    for f in range(frames):
        pc_list.append(np.dot(pca_evecs_sorted[:, p], r_vec[:, f]))
    pc_np = np.array(pc_list)
    mean_pc = np.mean(pc_np)
    dx = pc_np - mean_pc
    dx0 = np.zeros(dx.shape)
    dx0.fill(dx[0])
    mean_pc2 = np.mean(np.power(dx, 2))
    autoc = dx[0] * dx0[0] / mean_pc2
    print(dx[0], dx0[0], mean_pc2, autoc)

# plt.plot(range(N), np.power(pca_evecs_sorted[:, pcOfInterest], 2))
# plt.xticks(range(N), r_names)
# plt.show()

# plt.plot(range(frames), pc_list)
# plt.show()
