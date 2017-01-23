# location/paths for test/example data

from pkg_resources import resource_filename as res

DCD = res(__name__, './short.dcd')
PSF = res(__name__, './rpn11_ubq_interface-ionized.psf')

del res
