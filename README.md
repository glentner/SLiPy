# [SLiPy](http://glentner.github.io/SLiPy)

#### A Spectroscopy and astrophysical Library for Python 3

This Python package is an expanding code base for doing computational
astronomy, particularly spectroscopy. It contains both a *Spectrum* class
for handling spectra as objects (with +, -, \*, /, etc... operations defined)
and a growing suite of analysis tools.

**Dependencies:**
Python 3.x,
[astropy](http://www.astropy.org),
[matplotlib](http://matplotlib.org),
[numpy](http://www.numpy.org),
[scipy](http://www.scipy.org)

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![GitHub license](http://img.shields.io/badge/license-GPLv3-blue.svg?style=flat)](http://www.gnu.org/copyleft/gpl.html)
[![downloads](https://img.shields.io/github/downloads/glentner/slipy/latest/total.svg)](http://github.com/glentner/slipy)

Quick note: the subpackage **astrolibpy** was not developed
by me. It was coded by Sergey Koposov (@segasai) at Cambridge (then at least).
I found it useful for performing velocity corrections on my spectroscopic
data. I've modified several modules such that it can be imported and used in
Python 3.x. See his README file.

## Contents

SLiPy is split into several components. The principle component is the
subpackage **SLiPy** itself, which contains all the relevant
functionality. Further, **Data** is a package I'm working on that will provide
an API for searching astronomical data archives in a simple way. The other two
subpackages **Framework** and **astrolibpy** are of utility to the project but
not necessarily intended for export. As stated previously, astrolibpy was not
developed by me, only modified. I'm not going to document it's usage here. Its
name is unfortunate for me as it is a bit over done with the convention I was
already using, but for consistency I will keep it as it was from the author.

The following modules are elevated to the package level and are available
to import:

| SLiPy/ | Functions/Classes |
|---------|-------------------|
|[**Fits**](#FitsLoc)|[Find](#FindLoc), [RFind](#RFindLoc), [GetData](#GetDataLoc), [Header](#HeaderLoc), [Search](#SearchLoc), [PositionSort](#PositionSortLoc), |
|[**DataType**](#DataTypeLoc)|[WaveVector](#WaveVectorLoc), [Spectrum](#SpectrumLoc), |
|[**Simbad**](#SimbadLoc)|[Position](#PositionLoc), [Distance](#DistanceLoc), [Sptype](#SptypeLoc), [IDList](#IDListLoc), |
|[**Correlate**](#CorrelateLoc)|[XCorr](#XCorrLoc), |
|[**Telluric**](#TelluricLoc)|[Correct](#CorrectLoc), |
|[**Velocity**](#VelocityLoc)|[HelioCorrect](#HelioCorrectLoc), [BaryCorrect](#BaryCorrectLoc), [IrafInput](#IrafInputLoc),  |
|[**Observatory**](#ObservatoryLoc)|[OHP](#OHPLoc), |
|[**Plot**](#PlotLoc)|[SPlot](#SPlotLoc), [Iterate](#IterateLoc), |
|[**Montage**](#MontageLoc)|[SolveGrid](#MSolveGridLoc), [Mosaic](#MMosaicLoc), [SubField](#MSubFieldLoc), [Field](#MFieldLoc), |

| SLiPy/Data | Functions/Classes |
|------------|-------------------|
|[**Elodie**](#ElodieLoc)|[Archive](#EArchiveLoc), [Script](#EScriptLoc), [Download](#EDownloadLoc), |

##Installation

To install SLiPy, there is no setup procedure. Simply download the package,
either by clicking on the download link for a *tar* or *zip* archive or by
cloning it - `git clone http://github.com/glentner/SLiPy`. Extract it's
contents to wherever you like in a directory (ostensibly named *slipy*, but
actually you can call this library whatever you want as well because all
the imports are *relative*). Then add the parent directory to your `PYTHONPATH`
if it isn't already. For example:

```
$ cd
$ git clone http://github.com/glentner/SLiPy
$ echo "export PYTHONPATH=$PYTHONPATH:~" >> ~/.bashrc
```

And your ready to go!

##Exceptions

SLiPy attempts to catch all foreseeable exceptions and re-throw them under a
common handle with a human readable message. There is a unique exception class
for every module derived from *Exception*. The naming convention is for a
module's exception to be named after the module with the addition of the word
*Error*. So the *Fits* module will throw a *FitsError*.

##Contribute

If you use SLiPy or have your own code related to spectroscopy or computing
for astronomy and think it would be a useful addition (or you find a
bug/mistake) I'm more than open to suggested contributions/additions.

##Author

Geoffrey Lentner, B.S.  
Graduate Research Assistant  
Department of Physics & Astronomy  
University of Louisville

Website: [glentner.github.io](http://glentner.github.io)

##Acknowledge SLiPy

If you have made use of SLiPy in your project/research, you can acknowledge
your use in the following ways:

**Publications**  
This research has made use of SLiPy, an open source spectroscopy and
astrophysics library for Python 3 (G. Lentner, 2015).

**Projects**  
If your code either makes use of, borrows from SLiPy, a good way to reference
this is with a [shield](http://shields.io) in your README file.

[![SLiPy](http://img.shields.io/badge/running-SLiPy-green.svg?style=flat)](http://glentner.github.io/SLiPy)

The above badge is generated using the following snippet
```
[![SLiPy](http://img.shields.io/badge/running-SLiPy-green.svg?style=flat)](http://glentner.github.io/SLiPy)
```
# Documentation

#<a name=FitsLoc></a>[Fits](SLiPy/Fits.py)

Manipulate FITS files. Import data into *Spectrum* objects. Filter results
by right ascension and declination. Grab header elements. Search for attributes
of the data such as distance, spectral type, etc.

<a name=FindLoc></a>
- **Find** (*toplevel* = './', *pattern* = '\*.fits'):

    Search for file paths below *toplevel* fitting *pattern*. Returns a list
    of string values.

<a name=RFindLoc></a>
- **RFind** (*toplevel* = './', *pattern* = '\*.fits'):

    Recursively search for file paths below *toplevel* fitting *pattern*.
    Returns a list of string values.

<a name=GetDataLoc></a>
- **GetData** ( \**files*, \*\**kwargs*):

	Import data from FITS `files`. Returns a list of *Spectrum* objects.

    |Options     | Defaults        | Descriptions
    |------------|-----------------|------------------------------------------|
    |*verbose*   | True            | display messages, progress               |
    |*toplevel*  | ''              | request import from directory *toplevel* |  
    |*pattern*   | '\*.fits'       | pattern matching with *toplevel*         |
    |*recursive* | False           | search recursively below *toplevel*      |
    |*wavecal*   | True            | fit wavelength vector to data            |
    |*crpix1*    | 'crpix1'        | reference pixel header keyword           |
    |*crval1*    | 'crval1'        | value at reference pixel                 |
    |*cdelt1*    | 'cdelt1'        | resolution (delta lambda)                |
    |*xunits*    | 'Angstrom'      | wavelength units (astropy.units)         |
    |*yunits*    | 'ergs cm-2 s-1' | units of the data                        |

<a name=HeaderLoc></a>
- **Header** ( *filename*, *keyword*, \*\**kwargs*):

    Retrieve *keyword* from FITS header in file *filename*.
    Return type depends on what is returned.

<a name=SearchLoc></a>
- **Search** ( \**files*, \*\**kwargs*):

    Extract object names from Fits *files* and use Simbad module
    to resolve the *attribute* (a required keyword argument)
    from the SIMBAD astronomical database. Currently available attributes
    are 'Position', 'Distance', 'Sptype', and 'IDList'. Returns a list of
    results (type depends on the values).

    | Options     | Defaults  | Descriptions                         |
    |-------------|-----------|--------------------------------------|
    | *verbose*   | True      | display messages, progress           |
    | *toplevel*  | None      | search under *toplevel* directory    |
    | *pattern*   | '\*.fits' | for files under *toplevel*           |
    | *recursive* | False     | search recusively under *toplevel*   |
    | *attribute* | None      | attribute to search for (no default) |

<a name=PositionSortLoc></a>
- **PositionSort** ( *center*, *radius*, \**files*, \*\**kwargs* ):

    Return a list of files from *files* that lie in a *radius* (in
    decimal degrees) from *center*, based on the *ra* (right ascension) and
    *dec* (declination).

    | Options   | Defaults | Descriptions                             |
    |-----------|----------|------------------------------------------|
    *ra*        | 'pos1'   | header element for right ascension       |
    *dec*       | 'pos2'   | header element for declination           |
    *obj*       | 'object' | header element for object id             |
    *raconvert* | True     | convert decimal hours to decimal degrees |
    *verbose*   | True     | display messages, progress               |
    *toplevel*  | None     | *toplevel* directory to look for files   |
    *recursive* | False    | search *toplevel*ly below *toplevel*     |
    *pattern*   |'\*.fits' | glob *pattern* for file search           |
    *useSimbad* | False    | use *Simbad* instead of header elements  |


#<a name=DataTypeLoc></a>[DataType](SLiPy/DataType.py)

Objects for representing astronomical data. Currently, this includes the
*Spectrum* class and it's helper function *WaveVector*.

<a name=WaveVectorLoc></a>
- **WaveVector** ( *rpix*, *rval*, *delt*, *npix* ):

    Construct numpy array of wavelength values where *rpix* is the reference
    pixel index, *rval* is the wavelength at reference pixel, *delt* is the
    resolutions (delta lambda), and *npix* is the length of desired array.

<a name=SpectrumLoc></a>
- class **Spectrum** ( *filename*, *wavelengths* = None, \*\**kwargs* ):

    The *Spectrum* class is a container for a *data* array and it's
    corresponding wavelength calibration, *wave*. These are accessed with
    .data and .wave, respectively. The data is read in from the file,
    *filename*. Alternatively, it can be initialized by a set of numpy
    arrays. To do this, *filename* can actually be the *data* array and
    if this is the case, *wavelengths* must be assigned an array of equal
    length containing the corresponding wavelength values.

    The following operations are defined:
    +, -, \*, /, +=, -=, \*=, /=. For each of these, if the second operand is
    a scalar, the operation is performed pixel-wise on the *data* array. If
    the other operand is also a *Spectrum* object, the RHS operand is
    *resampled* onto the same pixel space as the LHS. The domain of the RHS
    *wave* array must be entirely contained by or equivalent to the LHS
    *wave* domain.

    The *wave* and *data* arrays are given units 'a la *astropy.units*;
    if the spectrum is initialized via numpy arrays, the units are only
    applied if there are none currently.

    | Options   | Defaults       | Descriptions                   |
    |-----------|----------------|--------------------------------|
    | *wavecal* | True           | fit wavelength vector to data  |
    | *crpix1*  | 'crpix1'       | reference pixel header keyword |
    | *crval1*  | 'crval1'       | value at reference pixel       |
    | *cdelt1*  | 'cdelt1'       | resolution (delta lambda)      |
    | *xunits*  | 'Angstrom'     | units for wavelength array     |
    | *yunits*  | 'erg cm-2 s-1' | units for data array           |

    Member functions:

    - *resample* ( \**args*, \*\**kwargs* ):

        If given a single argument, it is taken to be a *Spectrum* object,
        and *self* is resampled onto the pixel space of the other spectrum.
        Otherwise, three arguments are expected. The first and second argument
        should define the lower and upper wavelength value of a domain,
        respectively. The third argument should be the number of elements
        (pixels) for the new domain. Think numpy.linspace().

        | Options | Defaults | Descriptions                         |
        |---------|----------|--------------------------------------|
        | *kind*  | 'linear' | passed to scipy.interpolate.interp1d |

    - *copy* ():

        Essentially a wrapper to *deepcopy()*. To say SpectrumA = SpectrumB
        implies that SpectrumA *is* SpectrumB. If you want to create a new
        spectrum *equal* to another, say SpectrumA = SpectrumB.copy()


#<a name=SimbadLoc></a>[Simbad](SLiPy/Simbad.py)

This module allows the user to query the SIMBAD astronomical database from
inside Python or shell commands/scripts. It's four current major functions
*Position*, *Distance*, *Sptype*, and *IDList* return real variables with
appropriate types ready for use.

As a shell script:

```
$ Simbad.py

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

<a name=PositionLoc></a>
- **Position** ( *identifier*, \*\**kwargs* ):

    Return the right ascension and declination in decimal degrees of
    *identifier* as a pair.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *parse*   | True        | parse return file from SIMBAD    |
    | *full*    | False       | return more detailed information |

    Example:
    ```python
    ra, dec = Simbad.Position('Sirius')
    ```

<a name=DistanceLoc></a>
- **Distance** ( *identifier*, \*\**kwargs* ):

    Return the distance in parsecs to *identifier*.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *parse*   | True        | parse return file from SIMBAD    |
    | *full*    | False       | return more detailed information |

    Example:
    ```python
    distance = Simbad.Distance('rigel kent')
    ```

<a name=SptypeLoc></a>
- **Sptype** ( *identifier*, \*\**kwargs* ):

    Return the spectral type of *identifier* as resolved by SIMBAD.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *parse*   | True        | parse return file from SIMBAD    |
    | *full*    | False       | return more detailed information |

    Example:
    ```python
    # returns 'B8IVn' (HD 87901 is Regulus)
    sptype = Simbad.Sptype('HD 87901')
    ```

<a name=IDListLoc></a>
- **IDList** ( *identifier*, \*\**kwargs* ):

    Return a list of alternate IDs for *identifier*.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *parse*   | True        | parse return file from SIMBAD    |
    | *full*    | False       | return more detailed information |

    Example:
    ```python
    other_names = Simbad.IDList('proxima centauri')
    ```

#<a name=CorrelateLoc></a>[Correlate](SLiPy/Correlate.py)

Correlation functions for astronomical data.

<a name=XCorrLoc></a>
- **XCorr** ( *spectrumA*, *spectrumB*, \*\**kwargs* ):

    The function returns an integer value representing the best shift within
    a *lag* based on the computed RMS of each configuration.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *lag*     | 25          | pixel range to shift over        |

#<a name=TelluricLoc></a>[Telluric](SLiPy/Telluric.py)

Removal of atmospheric absorption lines in spectra.

<a name=CorrectLoc></a>
- **Correct** ( *spectrum*, \**calibration*, \*\**kwargs* ):

	Perform a telluric correction on *spectrum* with one or more
	*calibration* spectra. If more than one calibration spectrum is
	provided, the one with the best fit after performing both a
	horizontal cross correlation and a vertical amplitude fit is used.
	The spectrum and all the calibration spectra must have the same
	number of pixels (elements). If a horizontal shift in the calibration
	spectra is appropriate, only the corresponding range of the spectrum
	is divided out!

    **Notice** Your spectra must be continuum normalized for this to work!

    | Options   | Defaults        | Descriptions                         |
    |-----------|-----------------|--------------------------------------|
    | *lag*     | 25              | pixel range to shift over            |
    | *range*   | (0.5, 2.0, 151) | numpy.linspace for amplitude fitting |

    ![example](Figures/HD192640.png)

    **Figure 1:** The above figure is an example of applying the
    **Telluric.Correct()** algorithm to a spectrum. In this case, six spectra of
    *Regulus* from the Elodie archive were used as calibration spectra.


#<a name=VelocityLoc></a>[Velocity](SLiPy/Velocity.py)

Radial velocity corrections for 1D spectra.

<a name=HelioCorrectLoc></a>
- **HelioCorrect** ( *observatory*, \**spectra*, \*\**kwargs* ):

    Perform heliocentric velocity corrects on *spectra* based on
    *observatory* parameters (*longitude*, *latitude*, *altitude*) and the
    member attributes, *ra* (right ascension), *dec* (declination), and *jd*
    (julian date) from the *spectra*.

    | Options    | Defaults        | Descriptions               |
    |------------|-----------------|----------------------------|
    | *verbose*  | False           | display messages, progress |

<a name=BaryCorrectLoc></a>
- **BaryCorrect** ( *observatory*, \**spectra*, \*\**kwargs* ):

    Perform barycentric velocity corrects on *spectra* based on
    *observatory* parameters (*longitude*, *latitude*, *altitude*) and the
    member attributes, *ra* (right ascension), *dec* (declination), and *jd*
    (julian date) from the *spectra*.

    | Options    | Defaults        | Descriptions               |
    |------------|-----------------|----------------------------|
    | *verbose*  | False           | display messages, progress |

<a name=IrafInputLoc></a>
- **IrafInput** ( \**files*, \*\**kwargs* ):

	Build an input file for IRAF's rvcorrect task.

	*files* should be a list of FITS file names to build the output table for.
	The user can optionally specify a *toplevel* directory to search for files
    under. If *outfile* is given, write results to the file.

    | Options     | Defaults  | Descriptions                          |
    |-------------|-----------|---------------------------------------|
    | *toplevel*  | None      | search *toplevel* directory for files |
    | *pattern*   | '\*.fits' | pattern matching under *toplevel*     |
    | *recursive* | False     | search recusively under *toplevel*    |
    | *outfile*   | None      | write lines to file named *outfile*   |


#<a name=ObservatoryLoc></a>[Observatory](SLiPy/Observatory.py)

Define observatory parameter similar to the IRAF task. All observatories
should follow the following pattern. The user can add as many as they like
to this module. I welcome suggestions.

<a name=OHPLoc></a>
```Python
class OHP(Observatory):
	"""
	The Observatoire de Haute-Provence, France.
	"""
	def __init__(self):
		self.name      = 'Observatoire de Haute-Provence'
		self.latitude  = 43.9308334 # degrees N
		self.longitude = 356.28667  # degrees W
		self.altitude  = 650        # meters
```

#<a name=PlotLoc></a>[Plot](SLiPy/Plot.py)

Convenient wrappers to matplotlib for plotting spectra. A *SPlot* simply
creates a handle to remember figure attributes, to quickly go from looking
at one spectra to another. One can also *overlay* spectra.

<a name=SPlotLoc></a>
- class **SPlot** ( *spectrum*, \*\**kwargs* ):

    Spectrum Plot - Create figure of *spectrum*.

    | Options  | Defaults   | Descriptions         |
    |----------|------------|----------------------|
    | *marker* | 'b-'       | marker for data      |
    | *label*  | 'spectrum' | name of object       |
    | *usetex* | False      | render with pdflatex |

    The following member functions call pyplot equivalent:  
    **xlim**, **ylim**, **xlabel**, **ylabel**, **title**, **legend**,
    **text**, **grid**, **close**.

    Here, when these function are called, the arguments are passed to
    matplotlib; however, these calls are remembered. So when you go to *draw*
    the figure again, you are back where you left off.

    - *draw*( ):

        Rebuild and render the figure.

    - *refresh*( ):

        Re-render (refresh) the figure, without clearing the axis.

    - *txtclear*( ):

        Clear all the calls to *text*( ) from the figure.

    - *xoffset*( *value* ):

        Switch either on or off (*value* = True | False) the horizontal offset.

    - *yoffset*( *value* ):

        Switch either on or off (*value* = True | False) the vertical offset.

    - *overlay*( \**splots* ):

        Given one or more *splots*, *overlay* the figures.

    - *restore*( ):

        Restore the figure to only the original spectrum.

    - *save*( *filename* ):

        Save the figure to a file called *filename*. The file format is derived
        from the extention of *filename*.

<a name=IterateLoc></a>
- **Iterate**( \**splots*, \*\**kwargs* ):

	Iterate thru *splots* to inspect data, the user marks spectra of
	interest. The function returns a list of *keepers*.

    | Options | Defaults  | Descriptions          |
    |---------|-----------|-----------------------|
    | *keep*  | 'name'    | alternative is 'plot' |


#<a name=MontageLoc></a>[Montage](SLiPy/Montage.py)

[Montage](http://montage.ipac.caltech.edu/) is a very powerful suite of C code
for creating image mosaics. This module is a wrapper to automate not only the
process for small projects, but to segment large fields into a grid and
mosaic each *site* before combining them into a *master* mosaic.

<a name=MSolveGridLoc></a>
- **SolveGrid** ( *sides*, *grid* ):

	Helper function for the Field and SubField classes. Both *sides* and *grid*
	need to be array-like and of length two. *sides* is the side length of the
	field in decimal degrees in right ascension and declination respectively.
	*grid* specifies the subdivision along these axis (e.g., (2,2) says 2x2).

	The user should mindful of their choices. If the side lengths cannot be
	subdivided into well-behaved (rational) segments, higher decimal places
	will be lost in the SubField.ArchiveList() task resulting in small
	gaps in the mosaic.

<a name=MMosaicLoc></a>
- **Mosaic** ( *resolution*, \**folders*, \*\**kwargs* ):

    Conduct standard build procedures for all *folders*, similar to the
    [m101 example](http://montage.ipac.caltech.edu/docs/m101tutorial.html).
    *resolution* is the number of pixels per degree for the output image.
    Note: *folders* should be absolute paths. Further, below each of these
    directories, there should already exist the standard folder structure
    ```
    folder/
      |--images/
      |    |-- <location of raw FITS images>
      |
      |--projected/
      |--differences/
      |--corrected/
      |--final/
    ```

    The mosaic will be deposited at *final/mosaic.fits*.

    | Options    | Defaults | Descriptions                             |
    |------------|----------|------------------------------------------|
    | *verbose*  | True     | display messages, progress               |
    | *bkmodel*  | True     | model and correct for background effects |


<a name=MSubFieldLoc></a>
- class **SubField** ( *center*, *sides*, *grid*, \*\**kwargs* ):

    Create a grid of *sites* each of which will be mosaiced separately and
    then combined. Each of *center*, *sides*, and *grid* should be array-like
    and of length two. *center* should be the very center location for the
    mosaic in right ascension and declination (both in decimal degrees),
    respectively. *sides* needs to give the side lengths of the desired
    mosaic in decimal degrees (width-RA, height-DEC). *grid* should be the
    grid division for the field (e.g., (2, 2) means 2x2 grid).

    There will be a directory created for each *site* and also another *master*
    directory. The final resulting mosaic will be deposited at
    *final/mosaic.fits* here.

    | Options   | Defaults | Descriptions                             |
    |-----------|----------|------------------------------------------|
    | *verbose* | True     | display messages, progress               |
    | *survey*  | 'DSS'    | DSS, SDSS, 2MASS                         |
    | *band*    | 'DSS2B'  | filter for *survey*, see *bands* dict    |
    | *pad*     | 0.0      | amount to add (degrees) around *sites*   |

    The available filter band for each survey are as follows

    | Survey  | Bands                                                   |
    |---------|---------------------------------------------------------|
    | 'DSS'   | 'J', 'H', 'K'                                           |
    | 'SDSS'  | 'U', 'G', 'R', 'I', 'Z'                                 |
    | '2MASS' | 'DSS1B', 'DSS1R', 'DSS2B', 'DSS2R', 'DSS2IR', 'Quick-V' |

    The user should execute the following available methods in this order:

    - *ArchiveList* ( \*\**kwargs* ):

        Run the `mArchiveList` command on the *site* grid. The only keyword
        argument is *verbose* which defaults to True.

    - *ArchiveExec* ( \*\**kwargs* ):

        Run `mArchiveExec` on each *site* in the SubField. The only keyword
        argument is *verbose* which defaults to True.

    - *Build* ( *resolution*, \*\**kwargs*):

        Run the build process for the *sites* in this SubField. See the
        Montage.Mosaic() function documentation.

    - *Collect* ( \*\**kwargs* ):

        Collect the mosaics from all *site* locations into a master *images*
        folder. The only keyword argument is *verbose* which defaults to True.

    - *Merge* ( *resolution*, \*\**kwargs* ):

        Merge all *site* mosaics into a single master SubField mosaic. We are
        now calling Montage.Mosaic() on each *site*.

    **Example:**
    ```python
    from slipy import Montage, Simbad

    mosaic = Montage.SubField(

        Simbad.Position('Pleiades'), # center of the Pleiades
        (4, 4),                      # side lengths of mosaic
        (2, 2)                       # grid pattern creates 2x2 degree `sites`
    )

    mosaic.ArchiveList()
    mosaic.ArchiveExec()
    mosaic.Build()
    mosaic.Collect()
    mosaic.Merge()

    # If you have APLpy
    import aplpy

    fig = aplpy.FITSFigure('master/final/mosaic.fits')
    fig.show_grayscale()

    ```

<a name=MFieldLoc></a>
- class **Field** (*center*, *sides*, *grid*, *subgrid*, \*\**kwargs* ):

    Large image mosaic manager for Montage. This class (in terms of its
    usage) is the same as the *SubField* class, except that here we managing
    the subfields. So all the member functions are the same name and purpose,
    but instead act to call that same function on each subfield. Here, all
    the constructor arguments are as before, with the additions of *subgrid*
    which is also to be array-like of length two. *grid* will be the first
    level division to find the centers and side lengths of all the subfields
    and *subgrid* will be the further sub-division passed down to the daughter
    subfields.

    | Options   | Defaults | Descriptions                             |
    |-----------|----------|------------------------------------------|
    | *verbose* | True     | display messages, progress               |
    | *survey*  | 'DSS'    | DSS, SDSS, 2MASS                         |
    | *band*    | 'DSS2B'  | filter for *survey*, see *bands* dict    |
    | *pad*     | 0.0      | amount to add (degrees) around *sites*   |

    The available filter band for each survey are as follows

    | Survey  | Bands                                                   |
    |---------|---------------------------------------------------------|
    | 'DSS'   | 'J', 'H', 'K'                                           |
    | 'SDSS'  | 'U', 'G', 'R', 'I', 'Z'                                 |
    | '2MASS' | 'DSS1B', 'DSS1R', 'DSS2B', 'DSS2R', 'DSS2IR', 'Quick-V' |

    All the member functions are the same name as in *SubField*, but now with
    the addition of a final step:

    - *Finalize* ( *resolution*, \*\**kwargs* ):

        Collect all SubField/master mosaics into a single folder and
        run Mosaic() on them for a single final image.

#<a name=ElodieLoc></a>[Elodie](Data/Elodie.py)

Methods for data retrieval from the Elodie Archive.

<a name=EArchiveLoc></a>
- class **Archive** ( \*\**kwargs* ):

    Import and parse ascii catalog of Elodie archive files. The complete
    archive is stored in the member dictionary, *data*. It's organized
    by unique target names. Each target has a list of pairs consisting of the
	name of the file and the signal to noise for that spectrum. The reduced
    archive by default contains only *HD*, *BD*, *HR*, *GC*, and *GJ* objects,
    choosing the file pertaining to the spectra with the highest signal-to-noise
    ratio available.

    | Options    | Defaults                   | Descriptions          |
    |------------|----------------------------|-----------------------|
    | *infile*   | archives/elodie.csv        | path to input file    |
    | *catalogs* | ['HD','BD','HR','GC','GJ'] | catalogs to keep      |

    **Example:**
    ```python
    from slipy.Data import Elodie

    archive = Elodie.Archive()

    'HD187642' in archive.files # returns True (Altair)
    'HD045348' in archive.files # returns False (Canopus)
    ```

<a name=EScriptLoc></a>
- **Script** ( *filename*, *pipeline* = '' ):

    Construct url script for Elodie archive given *filename* and optionally
    *pipeline* instructions (e.g., '&z=wrs|fca[1,nor]').

<a name=EDownloadLoc></a>
- **Download** ( \**files*, \*\**kwargs* ):

    Download *files* from Elodie archive via url scripts. The spectra can be
    further reduced via Elodie's pipeline with the following options.

    | Options     | Defaults        | Descriptions                         |
    |-------------|-----------------|--------------------------------------|
    | *verbose*   | True            | display messages, progress           |
    | *resample*  | (min, max, res) | resample spectra (no default)        |
    | *normalize* | True            | continuum normalization              |
    | *outpath*   | './'            | directory for downloaded *files*     |
    | *names*     | []              | alternative output names for *files* |

    **Example:**
    ```python
    # all files in the archive for Altair (file name is first element in pair)
    files = [ x[0] for x in archive.data['HD187642'] ]

    # download files to current directory, resample spectrum to wavelengths
    # between 5850 and 5950 with a resolution of 0.01 Angstroms per pixel
    Elodie.Download( *files, resample=(5850, 5950, 0.01) )
    ```
