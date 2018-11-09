import dask.distributed
lc = dask.distributed.LocalCluster(n_workers=8, processes=True, memory_limit=200000000000)
client = dask.distributed.Client(lc)

from PyContact.core.ContactManager import ContactManager
from PyContact.core.ContactAnalyzer import ContactAnalyzer
# topo = "/mnt/drive/workspace/pycontactData/str_biot/nowater.psf"
# trajs = ["/mnt/drive/workspace/pycontactData/str_biot/trajectory-1.dcd"]
topo = "PyContact/exampleData/rpn11_ubq.psf"
trajs = ["PyContact/exampleData/rpn11_ubq.dcd"]
man = ContactManager(topo, trajs, 5.0, 2.5, 120, "segid RN11", "segid UBQ")
man.readTrajectories(nproc=8, use_pmda=True)

def summarize(man):
    import pandas as pd
    import numpy as np
    cs = man.atomicContactTrajectories[0].contacts
    cs_r = np.concatenate(cs)
    df = pd.DataFrame(cs_r, columns=["frame", "idx1", "idx2", "weight", "hbonds"])

    df["idx1"] = np.asarray(df["idx1"], dtype=np.int)
    df["idx2"] = np.asarray(df["idx2"], dtype=np.int)
    df["frame"] = np.asarray(df["frame"], dtype=np.int)
    df["hbonds"] = np.asarray(df["hbonds"], dtype=np.int)

    df["resname1"] = man.atomicContactTrajectories[0].resname_array[df["idx1"]]
    df["resname2"] = man.atomicContactTrajectories[0].resname_array[df["idx2"]]

    df["resid1"] = man.atomicContactTrajectories[0].resid_array[df["idx1"]]
    df["resid2"] = man.atomicContactTrajectories[0].resid_array[df["idx2"]]
    # df["accumulatePattern"] = df["resname1"] + df["resid1"].astype(str) + "-" + df["resname2"] + df["resid2"].astype(str)
    df2 = df.groupby(['frame', 'resname1', 'resid1', 'resname2', 'resid2'])['weight'].agg('sum').reset_index()
    return df2

summary = summarize(man)
