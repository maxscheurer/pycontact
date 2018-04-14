
class ContactManager(object):
    """
    Main manager for all contact data and trajectories


    """
    def __init__(self):
        trajectories = []



class ContactTrajectory:
    def __init__(self, keys, contactScores, bbScores, scScores,
                 hbonds):
        """
        Data container for trajectory contact data

        Parameters
        ----------
        keys: np.ndarray
            array of contact keys
        contactScores: np.ndarray
            accumulated contact scores per frame, shape(contact,frame)
        bbScores: tuple of 2 np.ndarrays
            backbone scores
        scScores: tuple of 2 np.ndarrays
            side chain scores
        hbonds: np.ndarray
            hydrogen bonds, shape(contact,frame)

        """
        self.keys = keys
        self.contactScores = contactScores
        self.bbScores1, self.bbScores2 = bbScores
        self.scScores1, self.scScores2 = scScores
        self.hbonds = hbonds
        self.numberOfFrames = len(self.keys)
