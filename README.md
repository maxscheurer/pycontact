# PyContact
Tools for contact analysis of biomolecules from Molecular Dynamics simulations.

## Installation
### Download
Download directly from github.

### Python 2.7
Unfortunately, pip2 does not provide binaries
for PyQt5, so it has to be [built from scratch](http://pyqt.sourceforge.net/Docs/PyQt5/installation.html).
Afterwards, run `python setup.py install`
to build and install PyContact.

### Python 3.5
Simply run `python setup.py install`.

**Dependency list**:
  * Qt5
  * PyQt5
  * numpy
  * matplotlib (built with qt5agg)
  * [MDAnalysis](http://www.mdanalysis.org)


* for visualization: latest version of [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)


## Screenshots
* Main Timeline Window
![timeline](screenshots/first_view_timeline.png?raw=true "Timeline view without filter")

* Filtered Timeline
![filtered](screenshots/filtered_contacts.png?raw=true "Filtered contacts")

* Preferences Panel
![pref](screenshots/pref_panel.png?raw=true "Preferences")

* Histogram Panel
![hist](screenshots/histogram_panel_new.png?raw=true "Histogram export")

## About
Authors: Maximilian Scheurer and Peter Rodenkirch

Please give a reference to this website if you use the software in research.
