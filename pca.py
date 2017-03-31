from __future__ import print_function
from PyContact.core.Scripting import (PyContactJob, JobConfig)
from PyContact.core.ContactFilters import ScoreFilter
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

jobAllMD = PyContactJob("./PyContact/exampleData/rpn11_ubq_interface-ionized.psf","./PyContact/exampleData/short.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,0,1,1,0], [0,0,0,1,1,0], "segid UBQ", "segid RN11"))

jobAllMD.runJob(8)
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
    for f in range(frames):
        r_vec[idx1, f] = c1.scoreArray[f]
        r_names[idx1] = c1.human_readable_title()
    idx1 += 1

pca = PCA(n_components=10)
components = pca.fit_transform(r_vec)
# print(components)
# plt.plot(components[:,0], components[:,1], 'o')
# plt.show()

model = TSNE(n_components=2, random_state=0)
results = model.fit_transform(components)
plt.plot(results[:,0], results[:,1], 'o')
plt.show()

