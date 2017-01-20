from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("wrap_cy",["wrap.pyx"],language="c++",
    extra_compile_args=["-std=c++11"],
    extra_link_args=["-std=c++11"]),
])
