""" pycontact setup
by Maximilian Scheurer, Peter Rodenkirch
"""
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


extensions = [Extension("PyContact.cy_modules.cy_gridsearch",
                        ["PyContact/cy_modules/cy_gridsearch.pyx"],
                        language="c++",
                        include_dirs=[".", "PyContact/cy_modules/src"],
                        extra_compile_args=["-std=c++0x"]), ]
setup(
    name='pycontact',
    version='1.0.5',
    description='PyContact',
    long_description='Tool for analysis of non-covalent interactions in MD trajectories',
    url='https://github.com/maxscheurer/pycontact',
    author='Maximilian Scheurer, Peter Rodenkirch',
    author_email='mscheurer@ks.uiuc.edu',
    license='GPLv3',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='computational biophysics simulation biology bioinformatics visualization protein biomolecules dna',
    package_dir={'PyContact': 'PyContact'},
    packages=find_packages(),
    python_requires=">=3.6",
    setup_requires=['cython'],
    install_requires=['numpy >= 1.16',
                      'matplotlib',
                      'mdanalysis >= 0.20.0',
                      'cython',
                      'seaborn',
                      'scipy',
                      'PyQt5'],
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(extensions),
    package_data={'PyContact': ['exampleData/defaultsession',
                                'exampleData/*.psf', 'exampleData/*.pdb',
                                'exampleData/*.dcd', 'exampleData/*.tpr',
                                'exampleData/*.xtc', 'gui/*.tcl',
                                'db/aa.db', 'cy_modules/*.pyx',
                                'cy_modules/src/*']},
    entry_points={
        'console_scripts': [
            'pycontact=PyContact.pycontact:main',
        ],
    },

)
