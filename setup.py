""" pycontact setup

by Maximilian Scheurer, Peter Rodenkirch
"""
from setuptools import setup, find_packages

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

    package_dir = {'pycontact': 'pycontact'},
    py_modules=['MainWindow'],
    packages=find_packages(),

    install_requires=['numpy','matplotlib','mdanalysis'],

    dependency_links = ['https://github.com/pyqt/python-qt5'],

    entry_points={
        'console_scripts': [
            'pycontact=pycontact:main',
        ],
    },

)
