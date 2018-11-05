from copy import deepcopy

import numpy as np

from pmda.parallel import ParallelAnalysisBase

from .ParallelAnalysis_core import loop_trajectory_grid


class ConvBond(object):
    """Python object of MDAnalysis bond for running jobs in parallel."""
    def __init__(self, bonds):
        super(ConvBond, self).__init__()
        self.types = []
        self.indices = []
        for b in bonds:
            self.indices.append(deepcopy(b.indices))
            self.types.append(deepcopy(b.type))

    def types(self):
        return self.types

    def to_indices(self):
        return self.indices


class AtomicContacts(ParallelAnalysisBase):
    """
        Calculate contacts
    """

    def __init__(self, universe,
                 sel1, sel2, config):
        """
        Parameters
        ----------
        universe: Universe
            MDA Universe
        sel1 : string
            Selection text 1
        sel2 : string
            Selection text 2
        config: list
            cutoff, hbondcutoff, hbondcutangle
        """
        self.config = config
        self.do_self_interaction = False

        if sel2 is "self":
            sel2 = sel1
            self.do_self_interaction = True

        atoms_sel1 = universe.select_atoms(sel1)
        atoms_sel2 = universe.select_atoms(sel2)
        super(AtomicContacts, self).__init__(universe,
                                             (atoms_sel1,
                                              atoms_sel2))
        atoms_all = universe.select_atoms("all")
        self.bonds = [ConvBond(a.bonds) for a in atoms_all.atoms]
        self.all_names = atoms_all.names
        self.all_resids = atoms_all.resids
        self.all_resnames = atoms_all.resnames
        self.all_segids = atoms_all.segids

    def _prepare(self):
        self.allContacts = None

    def _single_frame(self, ts, atomgroups):
        sel1, sel2 = atomgroups

        res = loop_trajectory_grid(sel1.positions,
                             sel2.positions,
                             sel1.indices,
                             sel2.indices,
                             self.config,
                             [self.bonds, self.all_names, self.all_resids, self.all_segids],
                             do_self_interaction=self.do_self_interaction, ts=ts.frame)

        return res

    def _conclude(self):
        import pandas as pd
        print("analysis done!")
        self.allContacts = np.hstack(self._results)
        self.result =  [self.allContacts, self.all_resnames,
                self.all_resids, self.all_names,
                self.all_segids, None]
