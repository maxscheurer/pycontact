""" pycontact setup
by Maximilian Scheurer, Peter Rodenkirch
"""
import os
import numpy
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# CUDA_FILES = ['PyContact/cy_modules/cy_sasa_cuda',
#               'PyContact/cy_modules/src/sasaCudaKernel',
#               ]

CUDA_FILES = ['PyContact/cy_modules/wrapper',
              'PyContact/cy_modules/src/manager',
              ]



def find_in_path(name, path):
    "Find a file in a search path"
    for dir in path.split(os.pathsep):
        binpath = os.path.join(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system
    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.
    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = os.path.join(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, or set $CUDAHOME')
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': os.path.join(home, 'include'),
                  'lib': os.path.join(home, 'lib')}
    for k, v in cudaconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))

    return cudaconfig


CUDA = locate_cuda()


extensions = [Extension("PyContact.cy_modules.cy_gridsearch",
                        ["PyContact/cy_modules/cy_gridsearch.pyx"], language="c++",
                        include_dirs=[".", "PyContact/cy_modules/src"]),
              Extension("PyContact.cy_modules.wrap_vmd",
                        ["PyContact/cy_modules/wrap_vmd.pyx"], language="c++",
                        include_dirs=[".", "PyContact/cy_modules/src"]),
              # Extension('PyContact.cy_modules.cy_sasa_cuda',
              #   sources=['PyContact/cy_modules/src/sasaCudaKernel.cu', 'PyContact/cy_modules/cy_sasa_cuda.pyx'],
              #   library_dirs=[CUDA['lib']],
              #   libraries=['cudart'],
              #   language='c++',
              #   runtime_library_dirs=[CUDA['lib']],
              #   extra_compile_args={'gcc': [],
              #                        'nvcc': ['-Wno-deprecated-gpu-targets', '-x=cu']},
              #   include_dirs = [CUDA['include'], 'PyContact/cy_modules/src']),
              Extension('PyContact.cy_modules.wrapper',
                sources=['PyContact/cy_modules/src/manager.cu', 'PyContact/cy_modules/wrapper.pyx'],
                library_dirs=[CUDA['lib']],
                libraries=['cudart'],
                language='c++',
                runtime_library_dirs=[CUDA['lib']],
                extra_compile_args={'gcc': [],
                                     'nvcc': ['-Wno-deprecated-gpu-targets', '-x=cu']},
                include_dirs = [CUDA['include'], 'PyContact/cy_modules/src']),
              ]


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.

    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""

    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        # if os.path.splitext(src)[1] == '.cu':
        if os.path.splitext(src)[0] in CUDA_FILES:
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # super(obj, src, ext, cc_args, extra_postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)

setup(
    name='pycontact',
    version='0.1.0b',
    description='pycontact - a tool for contact analysis of biomolecules from MD trajectories',
    long_description='',
    url='https://github.com/maxscheurer/pycontact',
    author='Maximilian Scheurer, Peter Rodenkirch',
    author_email='mscheurer@ks.uiuc.edu',
    license='GPLv3',

# https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
    ],

    keywords='computational biophysics simulation biology bioinformatics visualization protein biomolecules dna',

    package_dir = {'PyContact': 'PyContact'},
    packages=find_packages(),

    install_requires = ['numpy','matplotlib','mdanalysis','sklearn'],
    ext_modules=cythonize(extensions),
    cmdclass = {'build_ext': custom_build_ext},

    package_data = {'PyContact': ['exampleData/defaultsession','core/testpar.prm','exampleData/*.psf','exampleData/*.pdb','exampleData/*.dcd','exampleData/*.tpr','exampleData/*.xtc','gui/*.tcl','db/aa.db']},
    entry_points={
        'console_scripts': [
            'pycontact=PyContact.pycontact:main',
        ],
    },

)
