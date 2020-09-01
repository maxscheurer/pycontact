""" pycontact setup
by Maximilian Scheurer, Peter Rodenkirch
"""
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


extensions = [Extension("PyContact.cy_modules.cy_gridsearch",
                        ["PyContact/cy_modules/cy_gridsearch.pyx"], language="c++",
                        include_dirs=[".", "PyContact/cy_modules/src"], extra_compile_args=["-std=c++0x"]), ]
setup(
    name='pycontact',
    version='1.0.4',
    description='PyContact',
    long_description='Tool for analysis of non-covalent interactions in MD trajectories',
    url='https://github.com/maxscheurer/pycontact',
    author='Maximilian Scheurer, Peter Rodenkirch',
    author_email='mscheurer@ks.uiuc.edu',
    license='GPLv3',

# https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='computational biophysics simulation biology bioinformatics visualization protein biomolecules dna',

    package_dir = {'PyContact': 'PyContact'},
    packages=find_packages(),

    setup_requires = ['cython'],
    install_requires = ['numpy','matplotlib','mdanalysis','cython','seaborn', 'scipy'],
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(extensions),

    package_data = {'PyContact': ['exampleData/defaultsession','exampleData/*.psf','exampleData/*.pdb','exampleData/*.dcd','exampleData/*.tpr','exampleData/*.xtc','gui/*.tcl','db/aa.db','cy_modules/*.pyx','cy_modules/src/*']},
    entry_points={
        'console_scripts': [
            'pycontact=PyContact.pycontact:main',
        ],
    },

)
