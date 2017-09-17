from __future__ import print_function
from PyContact.core.Scripting import (PyContactJob, JobConfig)
from PyContact.core.ContactFilters import ScoreFilter
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd

n_comp = 2
jobAllMD = PyContactJob("./PyContact/exampleData/rpn11_ubq.psf","./PyContact/exampleData/rpn11_ubq.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid RN11", "segid UBQ"))
#
#
# jobAllMD = PyContactJob("/mnt/drive/workspace/fgf_dock_domain/test.pdb","/mnt/drive/workspace/fgf_dock_domain/fgf2_a1_domainrosetta_scoring_min_1500.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,0,0,1], [0,0,1,1,0], "segid A", "segid B"))


# jobAllMD = PyContactJob("/mnt/drive/workspace/fgf_atp1a1_human_docking/dock_fgf2-cluster_rep_1.pdb","/mnt/drive/workspace/fgf_atp1a1_human_docking/atp1a1_fgf2_complexrosetta_scoring_min_5000.dcd", "test", JobConfig(5.0, 2.5, 120, [0,0,1,1,0], [0,0,1,1,0], "segid A", "segid B"))

jobAllMD.runJob(8)
# jobAllMD.writeSessionToFile()
contacts = jobAllMD.analyzer.finalAccumulatedContacts

filter = ScoreFilter("score", "greater", 0.1, u"Median")
contacts = filter.filterContacts(contacts)

N = len(contacts)
frames = len(jobAllMD.analyzer.contactResults)
print("frms:", frames, "N_contacts: ", N)

raw_data = {}
for c1 in contacts:
    current_list = np.array(c1.scoreArray)
    # current_list[current_list > 0] = 1
    raw_data[c1.human_readable_title()] = current_list
df = pd.DataFrame(raw_data)

df = df.transpose()
print(df.head())

pca = PCA(n_components=n_comp)
pca.fit(df)
pca_2d = pca.transform(df)
df_pca = pd.DataFrame(pca_2d)
df_pca.index = df.index
df_pca.columns = ['PC1','PC2']
print(pca.explained_variance_ratio_)
ax = df_pca.plot(kind='scatter', x="PC2", y="PC1", figsize=(16,8))
for f, c in enumerate(df.index):
    ax.annotate(c, (df_pca.iloc[f].PC2, df_pca.iloc[f].PC1))
plt.show()
