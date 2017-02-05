from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("wrap_vmd",["wrap_vmd.pyx"],language="c++"
    ,include_dirs=['/usr/local/include']
),
])
