""" pycontact setup
by Maximilian Scheurer, Peter Rodenkirch
"""
from setuptools import setup, find_packages

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
        'Programming Language :: Python :: 3.7',
    ],
    keywords='computational biophysics simulation biology bioinformatics visualization protein biomolecules dna',
    package_dir={'PyContact': 'PyContact'},
    packages=find_packages(),
    python_requires=">=3.7",
    setup_requires=['cython'],
    install_requires=['numpy >= 1.16',
                      'matplotlib',
                      'mdanalysis >= 2.0.0',
                      'seaborn',
                      'scipy',
                      'PyQt5'],
    package_data={'PyContact': ['exampleData/defaultsession',
                                'exampleData/*.psf', 'exampleData/*.pdb',
                                'exampleData/*.dcd', 'exampleData/*.tpr',
                                'exampleData/*.xtc',
                                ]},
    entry_points={
        'console_scripts': [
            'pycontact=PyContact.pycontact:main',
        ],
    },

)
