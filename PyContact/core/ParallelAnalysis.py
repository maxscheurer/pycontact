import numpy as np

from pmda.parallel import ParallelAnalysisBase

from .ParallelAnalysis_core import loop_trajectory_grid


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
        self.universe = universe
        self.sel1 = sel1
        self.sel2 = sel2
        self.config = config
        self.do_self_interaction = False

        if self.sel2 is "self":
            self.sel2 = self.sel1
            self.do_self_interaction = True

        self.atoms_sel1 = self.universe.select_atoms(self.sel1)
        self.atoms_sel2 = self.universe.select_atoms(self.sel2)
        super(AtomicContacts, self).__init__(self.universe,
                                             (self.atoms_sel1,
                                              self.atoms_sel2))

    def _prepare(self):
        self.atoms_all = self.universe.select_atoms("all")
        self.atoms_bb = self.universe.select_atoms("backbone")
        self.allContacts = None

    def _single_frame(self, ts, atomgroups):
        sel1, sel2 = atomgroups

        loop_trajectory_grid(sel1.positions,
                             sel2.positions,
                             sel1.indices,
                             sel2.indices,
                             self.config,
                             [self.atoms_all.bonds,
                              self.atoms_all.names,
                              self.atoms_all.resid],
                             do_self_interaction=self.do_self_interaction)

        return None

    def _conclude(self):
        self.timeseries = np.hstack(self._results)
