# PyContact
[![Build Status](https://travis-ci.com/maxscheurer/pycontact.svg?token=Xyntx2ELmeydq8pgqs8t&branch=develop)](https://travis-ci.com/maxscheurer/pycontact)

**PyContact** is a GUI-based tool for rapid, versatile and customizable analysis of non-covalent interactions in molecular dynamics (MD) trajectories.

## Installation
### Download
Download directly from github.

### Python 2.7
Unfortunately, pip2 does not provide binaries
for PyQt5, so it has to be [built from scratch](http://pyqt.sourceforge.net/Docs/PyQt5/installation.html).
Afterwards, run `python setup.py install`
to build and install PyContact.

### Python 3.5
MDAnalysis does currently NOT support DCD files in Python 3 and is not well tested.


**Dependency list**:
  * Qt5
  * PyQt5
  * numpy
  * matplotlib (built with qt5agg)
  * [MDAnalysis](http://www.mdanalysis.org)


* for visualization: latest version of [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)


## Scripting
Analysis can also be performed by using PyContact as a plain python package. In this way,
you automate analysis (reading trajectories and running contact score accumulation) from your python script and
later on visualize it in the GUI.

## Citation
If you use PyContact in a scientific publication, please cite:
(placeholder)

## About
Main Developers: Maximilian Scheurer and Peter Rodenkirch
