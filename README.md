# AstroPython

#### A Computational Astronomy Package for Python 3

This is a Python package containing modules I've developed to speed 
up my work flow. It mostly just contains functions for doing spectroscopy. 
It will continue to be developed. I'm sharing it in the hope that others 
may find it helpful. 

**Dependencies:** Python 3.x, astropy, matplotlib, numpy, scipy

Quick note: the subpackage **astrolibpy** was not developed
by me. It was coded by Sergey Koposov (@segasai) at Cambridge (then at least). 
I found it useful for performing velocity corrections on my spectroscopic 
data. I've modified several modules such that it can be imported and used in 
Python 3.x. See his README file.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Modules:

AstroPython is split into several components. The principle component is the
subpackage *AstroPython* itself, which contains all the relavent functionality.
The other two subpackages *Framework* and *astrolibpy* are of utility to the
project but necessarily intended for export. As stated previously, *astrolibpy* 
was not developed by me, only modified. I'm not going to document it's usage 
here. Its name is unfortunate for me as it is a bit over done with the 
convention I was already using, but for consistency I will keep it as it was 
from the author.

The following modules are elevated to the package level and are available
to import whole:

- Fits 
- Simbad
- Correlate
- Telluric 
- Velocity
- Montage
- Observatory


Documentation on the specific tools available here is forthcoming. In the 
interim, most are self-documenting.
