# location/paths for test/example data

from pkg_resources import resource_filename as res

DCD = res(__name__, './rpn11_ubq.dcd')
PSF = res(__name__, './rpn11_ubq.psf')
XTC = res(__name__, './md_noPBC.xtc')
TPR = res(__name__, './md.tpr')
DEFAULTSESSION = res(__name__, './defaultsession')

del res
