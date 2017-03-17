class Configuration(object):
    """Sets the current configuration."""
    def __init__(self, psf, dcd, cutoff, hbondcutoff, hbondcutangle, sel1text, sel2text):
        super(Configuration, self).__init__()
        self.psf = psf
        self.dcd = dcd
        self.cutoff = cutoff
        self.hbondcutoff = hbondcutoff
        self.hbondcutangle = hbondcutangle
        self.sel1text = sel1text
        self.sel2text = sel2text
