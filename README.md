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

* [**Fits**](#Fits) [Find, RFind, GetData, ]

* [**DataType**](#DataType) [WaveVector, Spectrum, ]

* [**Simbad**](#Simbad) [Query, Position, Distance, Sptype, IDList, ]

* [**Correlate**](#Correlate) [XCorr, ]

* [**Telluric**](#Telluric) [Correct, ]

* [**Velocity**](#Velocity) [HelioCorrect, ]

* [**Montage**](#Montage)

* [**Observatory**](#Observatory) [OHP, ]


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

#<a name=Fits></a>[Fits](AstroPython/Fits.py)

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
#<a name=DataType></a>[DataType](AstroPython/DataType.py)

Objects for representing astronomical data.

```Python
def WaveVector( rpix, rval, delt, npix ):
	"""
	Construct numpy array of wavelength values based on:

		`rpix` : reference pixel index
		`rval` : wavelength at reference pixel
		`delt` : resolutions (delta lambda)
		`npix` : length of desired array
	"""
```
```Python
class Spectrum:
	"""
	Spectrum objects consist of a `data` vector, and optionally a  
	`wavelength` vector (accessed with .data and .wave respectively).
    (+, -, *, /, +=, -=, *=, /=) are overloaded. The LHS spectrum is the
    reference and the RHS spectrum is resampled onto the wavelength space
    of the LHS spectrum before applying operations pixel-wise. Scalar
    operations applied to all pixels.
	"""
	def __init__(self, argument, **kwargs ):
		"""
		Hold spectrum `data` from file name `argument`. Alternatively,
		construct spectra with a one-dimensional numpy.ndarray as `argument`.
		`wavecal` is assumed to be true for file input.

		kwargs = {
    		'wavecal': True     , # fit wavelength vector to data
    		'crpix1' : 'crpix1' , # reference pixel header keyword
    		'crval1' : 'crval1' , # value at reference pixel
    		'cdelt1' : 'cdelt1' , # resolution (delta lambda)
		}
    """
```

#<a name=Simbad></a>[Simbad](AstroPython/Simbad.py)

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

The following objects/functions are available:

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
```Python
def Position( identifier, **kwargs ):
	"""
	Position( identifier, **kwargs ):

	Handle to the Query class with criteria='%C00(d;C)'. Return right
    ascension and declination in decimal degrees of `identifier`.

    Example:

    ra, dec = Position('Sirius')
	"""
```
```Python
def Distance( identifier, **kwargs ):
	"""
	Distance( identifier, **kwargs ):

	Handle to the Query class with criteria='%PLX' Return the distance
    in parsecs to `identifier`.

    Example:

    d = Distance('rigel kent')
	"""
```
```Python
def Sptype(identifier, **kwargs):
	"""
	Sptype( identifier, **kwargs ):

	Handle to the Query class with criteria='%SP'. Return the
    spectral type as resolved by SIMBAD.

    Example:

    sptype = Sptype('HD 87901') # returns 'B8IVn' (HD 87901 is Regulus)
	"""
```
```Python
def IDList(identifier, **kwargs):
	"""
	IDList(identifier, **kwargs):

	Handle to the Query class with criteria='%IDLIST'.
	With `parse` = True, return a list of alternate IDs for
	the `identifier` provided.
	"""

    Example:

    other_names = IDList('proxima centauri')
```

#<a name=Correlate></a>[Correlate](AstroPython/Correlate.py)

Module of correlation functions for astronomical data.

```Python
def Xcorr( spectrumA, spectrumB, **kwargs ):
	"""
	Cross correlate two spectra of equal pixel length. The function returns
	an integer value representing the best shift within a `lag` based on
	the computed RMS of each configuration.
	"""
```

#<a name=Telluric></a>Telluric

Removal of atmospheric adsorption lines in spectra.

```Python
def Correct(spectrum, *calibration, **kwargs):
	"""
	Perform a telluric correction on `spectrum` with one or more
	`calibration` spectra. If more than one calibration spectrum is
	provided, the one with the best fit after performing both a
	horizontal cross correlation and a vertical amplitude fit is used.
	The spectrum and all the calibration spectra must have the same
	number of pixels (elements). If a horizontal shift in the calibration
	spectra is appropriate, only the corresponding range of the spectrum
	is divided out!

	kwargs = {
			lag  : 25            , # pixel shift limit for XCorr()
			range:(0.5, 2.0, 151), # numpy.linspace for amplitude trials
		}
	"""
```

#<a name=Velocity></a>[Velocity](AstroPython/Velocity.py)

#<a name=Observatory></a>[Observatory](AstroPython/Observatory.py)

#<a name=Montage></a>[Montage](AstroPython/Montage.py)
