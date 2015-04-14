# AstroPython

#### A Computational Astronomy Package for Python 3

This is a Python package containing modules I've developed to speed
up my work flow. It mostly just contains functions for doing spectroscopy.
It will continue to be developed. I'm sharing it in the hope that others
may find it helpful.

**Dependencies:** Python 3.x, astropy, matplotlib, numpy, scipy

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

Quick note: the subpackage **astrolibpy** was not developed
by me. It was coded by Sergey Koposov (@segasai) at Cambridge (then at least).
I found it useful for performing velocity corrections on my spectroscopic
data. I've modified several modules such that it can be imported and used in
Python 3.x. See his README file.

## Modules

AstroPython is split into several components. The principle component is the
subpackage **AstroPython** itself, which contains all the relavent
functionality. Further, **Data** is a package I'm working on that will provide
an API for searching astronomical data archives in a simple way. The other two
subpackages **Framework** and **astrolibpy** are of utility to the project but
not necessarily intended for export. As stated previously, astrolibpy was not
developed by me, only modified. I'm not going to document it's usage here. Its
name is unfortunate for me as it is a bit over done with the convention I was
already using, but for consistency I will keep it as it was from the author.

The following modules are elevated to the package level and are available
to import:

##[Fits](AstroPython/Fits.py)

Import data from, handle, and manipulate FITS format files.

```Python
def Find(toplevel = './', pattern = '*.fits'):
    """
    Search for file paths below `toplevel` fitting `pattern`. Returns a list
    of string values.
    """
```

```Python
def RFind(toplevel = './', pattern = '*.fits'):
    """
    Recursively search for paths below `toplevel` fitting `pattern`. Returns
    a list of string values.
    """
```

```Python
def GetData( *files, **kwargs ):
	"""
	Import data from FITS `files`. Returns a list of Spectrum objects.

	kwargs = {
			verbose   : True    , # display messages, progress
			toplevel  : ''      , # request import from directory `toplevel`
			pattern   : '*.fits', # pattern matching with `toplevel`
			recursive : False   , # search recursively below `toplevel`
			wavecal   : True    , # fit wavelength vector to data
			crpix1    : 'crpix1', # reference pixel header keyword
			crval1    : 'crval1', # value at reference pixel
			cdelt1    : 'cdelt1', # resolution (delta lambda)
		}
	"""
```

```Python
def Header( filename, keyword, **kwargs ):
	"""
	Retrieve `keyword` from FITS header in `filename`. Return type depends on
    the value.
	"""
```

```Python
def Search( *files, **kwargs ):
	"""
	Extract object names from Fits `files` and use Simbad module
	to resolve the `attribute` (a required keyword argument)
	from the SIMBAD astronomical database. Currently available attributes
    are 'Position', 'Distance', 'Sptype', and 'IDList'. Returns a list of
    results (type depends on the values).

	kwargs = {
			verbose   : True    , # display messages, progress
			toplevel  : ''      , # search under `toplevel` directory
			pattern   : '*.fits', # for files under `toplevel`
			recursive : False   , # search recusively under `toplevel`
			attribute : ''      , # attribute to search for (no default)
		}
	"""
```

```Python
def PositionSort( center, radius, *files, **kwargs ):
    """
    Return a list of files from `files` that lie in a `radius` (in degrees)
    from `center`, based on the `ra` (right ascension) and `dec` (declination).

    kwargs = {
            'ra'       : 'pos1'  , # header element for right ascension
            'dec'      : 'pos2'  , # header element for declination
            'obj'      : 'object', # header element for object id
            'raconvert': True    , # convert decimal hours to decimal degrees
            'verbose'  : True    , # display messages, progress
            'toplevel' : ''      , # `toplevel` directory to look for files in
            'recursive': False   , # search `recursive`ly below `toplevel`
            'pattern'  : '*.fits', # glob `pattern` for file search
            'useSimbad': False     # use Simbad instead of header elements
        }
    """

```

##[Simbad](AstroPython/Simbad.py)

This module allows the user to query the SIMBAD astronomical database from
inside Python or shell commands/scripts.

As a shell script:

```
usage: Simbad.py @Attribute <identifier> [**kwargs]

The 'Attribute' points to a function within this module and indicates
what is to be run. Execute 'Simbad.py @Attribute help' for usage details of
a specific function. Currently available attributes are: `Position`,
`Distance`, `Sptype` and `IDList`.

The identifier names can be anything recognized by SIMBAD (e.g., Regulus,
"alpha leo", "HD 121475", "del cyg", etc ...) if the name is two parts make
sure to use quotation to enclose it.

The **kwargs is the conventional reference to Python keyword arguments.
These should be specific to the 'Attribute' being pointed to.
```

The following objects, functions are available:

```Python
class Query:
	"""
	Query( identifier, criteria, **kwargs ):

	Class for querying the SIMBAD astronomical database for 'citeria'
	of 'identifier'. This object is not intended to be used directly; it
    is an abstraction and is used by the other functions which should be
    called by the user.

	kwargs = {
		'parse' : True,  # extract relevant data from SIMBAD return file
		'dtype' : float, # output datatype
	}
	"""
```

###[Correlate](AstroPython/Correlate.py)


###[Telluric](AstroPython/Telluric.py)


###[Velocity](AstroPython/Velocity.py)


###[Montage](AstroPython/Montage.py)


###[Observatory](AstroPython/Fits.py)


Documentation on the specific tools available here is forthcoming. In the
interim, most are self-documenting.
