# [SLiPy](http://glentner.github.io/SLiPy)

#### A Spectroscopy and astrophysics Library for Python 3

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
|[**Observatory**](#ObservatoryLoc)|[...](#OHPLoc), |
|[**Plot**](#PlotLoc)|[SPlot](#SPlotLoc), [Iterate](#IterateLoc), |
|[**Profile**](#ProfileLoc)|[Select](#SelectLoc), [Fit](#FitLoc), [Extract](#ExtractLoc),   |
|[**Montage**](#MontageLoc)|[SolveGrid](#MSolveGridLoc), [Mosaic](#MMosaicLoc), [SubField](#MSubFieldLoc), [Field](#MFieldLoc), |

<br>

| SLiPy/Data | Functions/Classes |
|------------|-------------------|
|[**Elodie**](#ElodieLoc)|[Archive](#EArchiveLoc), [Script](#EScriptLoc), [Download](#EDownloadLoc), |

## Installation

To install SLiPy, there is no setup procedure. Simply download the package,
either by clicking on the download link for a *tar* or *zip* archive or by
cloning it.

Extract it's contents to wherever you like in a directory (ostensibly named
*slipy*, but actually you can call this library whatever you want as well
because all the imports are *relative*). Then add the parent directory to your
*PYTHONPATH* if it isn't already. For example:

```
$ cd
$ git clone http://github.com/glentner/slipy
$ echo "export PYTHONPATH=$PYTHONPATH:~" >> ~/.bashrc
```

And your ready to go!

## Exceptions

SLiPy attempts to catch all foreseeable exceptions and re-throw them under a
common handle with a human readable message. There is a unique exception class
for every module derived from *Exception*. The naming convention is for a
module's exception to be named after the module with the addition of the word
*Error*. So the *Fits* module will throw a *FitsError*.

## Contribute

If you use SLiPy or have your own code related to spectroscopy or computing
for astronomy and think it would be a useful addition (or you find a
bug/mistake) I'm more than open to suggested contributions/additions.

## Author(s)

[Geoffrey Lentner](http://glentner.github.io), B.S.  
Graduate Research Assistant  
Department of Physics & Astronomy  
University of Louisville

Significant intellectual contributions have been made by my thesis advisor,
specifically in terms of the science behind much of this package.

[James Lauroesch](http://www.physics.louisville.edu/jtl/), Ph.D.  
Associate Professor  
Department of Physics & Astronomy  
University of Louisville


##Acknowledge SLiPy

If you have made use of SLiPy in your project/research, you can acknowledge
your use in the following ways:

**Publications**  
This research has made use of SLiPy, an open source spectroscopy and
astrophysics library for Python 3 (G. Lentner, 2015).

**Projects**  
If your code either makes use of or borrows from SLiPy, a good way to reference
this is with a [shield](http://shields.io) in your README file.

[![SLiPy](http://img.shields.io/badge/running-SLiPy-green.svg?style=flat)](http://glentner.github.io/SLiPy)

The above badge is generated using the following snippet
```
[![SLiPy](http://img.shields.io/badge/running-SLiPy-green.svg?style=flat)](http://glentner.github.io/SLiPy)
```

# Documentation

<br>
##<a name=FitsLoc></a>[Fits](SLiPy/Fits.py)

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

	Import data from FITS *files*. Returns a list of *Spectrum* objects.

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
    *recursive* | False    | search below *toplevel* recursively      |
    *pattern*   |'\*.fits' | glob *pattern* for file search           |
    *useSimbad* | False    | use *Simbad* instead of header elements  |


<br>
##<a name=DataTypeLoc></a>[DataType](SLiPy/DataType.py)

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

    <br>
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


<br>
##<a name=SimbadLoc></a>[Simbad](SLiPy/Simbad.py)

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

    <br>
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

    <br>
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

    <br>
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

    <br>
    Example:
    ```python
    other_names = Simbad.IDList('proxima centauri')
    ```


<br>
##<a name=CorrelateLoc></a>[Correlate](SLiPy/Correlate.py)

Correlation functions for astronomical data.

<a name=XCorrLoc></a>
- **XCorr** ( *spectrumA*, *spectrumB*, \*\**kwargs* ):

    The function returns an integer value representing the best shift within
    a *lag* based on the computed RMS of each configuration.

    | Options   | Defaults    | Descriptions                     |
    |-----------|-------------|----------------------------------|
    | *lag*     | 25          | pixel range to shift over        |


<br>
##<a name=TelluricLoc></a>[Telluric](SLiPy/Telluric.py)

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


<br>
##<a name=VelocityLoc></a>[Velocity](SLiPy/Velocity.py)

Radial velocity corrections for 1D spectra.

<a name=HelioCorrectLoc></a>
- **HelioCorrect** ( *observatory*, \**spectra*, \*\**kwargs* ):

    Perform heliocentric velocity corrects on *spectra* based on
    *observatory* parameters (*longitude*, *latitude*, *altitude*) and the
    member attributes, *ra* (right ascension), *dec* (declination), and *jd*
    (julian date) from the *spectra*. These should all have units.

    | Options    | Defaults        | Descriptions               |
    |------------|-----------------|----------------------------|
    | *verbose*  | False           | display messages, progress |

<a name=BaryCorrectLoc></a>
- **BaryCorrect** ( *observatory*, \**spectra*, \*\**kwargs* ):

    Perform barycentric velocity corrects on *spectra* based on
    *observatory* parameters (*longitude*, *latitude*, *altitude*) and the
    member attributes, *ra* (right ascension), *dec* (declination), and *jd*
    (julian date) from the *spectra*. These should all have units.

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


<br>
##<a name=ObservatoryLoc></a>[Observatory](SLiPy/Observatory.py)

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
        self.longitude = 356.28667  * u.degree # West
        self.latitude  = 43.9308334 * u.degree # North
        self.altitude  = 650        * u.meter
        self.timezone  = 1          * u.hourangle
```

There are currently 69 defined observatories:

| Class ID       | Observatory Name                                         |
|----------------|----------------------------------------------------------|
|            OHP |              The Observatoire de Haute-Provence, France. |
|           KPNO |                           Kitt Peak National Observatory |
|           WIYN |                                         WIYN Observatory |
|           CTIO |                   Cerro Tololo Interamerican Observatory |
|        LASILLA |                 European Southern Observatory: La Silla. |
|        PARANAL |                   European Southern Observatory: Paranal |
|           LICK |                                         Lick Observatory |
|           MMTO |                                          MMT Observatory |
|           CFHT |                           Canada-France-Hawaii Telescope |
|        LAPALMA |                        Roque de los Muchachos, La Palma. |
|            MSO |                                  Mt. Stromlo Observatory |
|            SSO |                                Siding Spring Observatory |
|            AAO |                             Anglo-Australian Observatory |
|       MCDONALD |                                     McDonald Observatory |
|            LCO |                                 Las Campanas Observatory |
|      MTBIGELOW |                  Catalina Observatory: 61 inch telescope |
|            DAO |                       Dominion Astrophysical Observatory |
|            SPM |     Observatorio Astronomico Nacional, San Pedro Martir. |
|           TONA |         Observatorio Astronomico Nacional, Tonantzintla. |
|        PALOMAR |                                       The Hale Telescope |
|            MDM |                       Michigan-Dartmouth-MIT Observatory |
|            NOV |                        National Observatory of Venezuela |
|            BMO |                              Black Moshannon Observatory |
|            BAO |                             Beijing XingLong Observatory |
|           KECK |                                   W. M. Keck Observatory |
|           EKAR |                               Mt. Ekar 182 cm. Telescope |
|            APO |                                 Apache Point Observatory |
|         LOWELL |                                       Lowell Observatory |
|            VBO |                                  Vainu Bappu Observatory |
|            IAO |                   Indian Astronomical Observatory, Hanle |
|           FLWO |                                      Whipple Observatory |
|          FLWO1 |                                      Whipple Observatory |
|            ORO |                                    Oak Ridge Observatory |
|            LNA |             Laboratorio Nacional de Astrofisica - Brazil |
|           SAAO |                   South African Astronomical Observatory |
|         CASLEO |              Complejo Astronomico El Leoncito, San Juan. |
|         BOSQUE |             Estacion Astrofisica Bosque Alegre, Cordoba. |
|         ROZHEN |     National Astronomical Observatory Rozhen - Bulgaria. |
|           IRTF |                         NASA Infrared Telescope Facility |
|          BGSUO |                    Bowling Green State Univ Observatory. |
|           DSAZ |     Deutsch-Spanisches Observatorium Calar Alto - Spain. |
|             CA |                                   Calar Alto Observatory |
|           HOLI |  Observatorium Hoher List (Universitaet Bonn) - Germany. |
|            LMO |                            Leander McCormick Observatory |
|            FMO |                                 Fan Mountain Observatory |
|         WHITIN |                     Whitin Observatory,Wellesley College |
|            OSN |                            Observatorio de Sierra Nevada |
|   GEMINI NORTH |                                 Gemini North Observatory |
|   GEMINI SOUTH |                                 Gemini South Observatory |
|        LASILLA |                 European Southern Observatory: La Silla. |
|        PARANAL |                  European Southern Observatory: Paranal. |
|         ESONTT |            European Southern Observatory, NTT, La Silla. |
|         ESO36M | European Southern Observatory, 3.6m Telescope, La Silla. |
|         ESOVLT |             European Southern Observatory, VLT, Paranal. |
|            SLN |                 SLN - Catania Astrophysical Observatory. |
|            EUO |                               Ege University Observatory |
|            TNO |                             Turkiye National Observatory |
|            TUG |                    TUBITAK National Observatory, Turkey. |
|            MGO |                                 Mount Graham Observatory |
|          ARIES | Aryabhatta Research Institute of Observational Sciences. |
|           OALP |                     Observatorio Astronomico de La Plata |
|           OLIN |                   Connecticut College - Olin Observatory |
|         BOYDEN |                                       Boyden Observatory |
|          LULIN |                                        Lulin Observatory |
|           SOAR |               Southern Astrophysical Research Telescope. |
|          BAKER |                                        Baker Observatory |
|            HET |           McDonald Observatory - Hobby-Eberly Telescope. |
|           JCDO |        Jack C. Davis Observatory, Western Nevada College |
|            LNO |                            Langkawi National Observatory |


<br>
##<a name=PlotLoc></a>[Plot](SLiPy/Plot.py)

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

    <br>
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
		
	- *markers*( \**args* ):
	
		Reassign the values for the `marker`s in the figure. The number
		of arguments must equal the number of spectra in the figure. This
		starts out as one, but will increase for ever call to SPlot.overlay().

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

    <br>
    **Example:**
    ```python
    from slipy import Fits, Plot

    fpath = '?' # toplevel directory name where your FITS files are
    files = Fits.Find(fpath)

    spectra = Fits.GetData( *files )

    figure = [

        Plot.SPlot( spectrum, label=Fits.Header(fname, 'object') )
        for spectrum, fname in zip(spectra, files)
    ]

    keepers = Plot.Iterate( *figure )
    # enter either 'y', 'n', or 'x' as prompted by the terminal
    ```


<br>
##<a name=ProfileLoc></a>[Profile](SLiPy/Profile.py)

Profile fitting tasks for spectra.

<a name=SelectLoc></a>
- **Select** ( *splot* ):

    Select points from the *splot*. This should be of type SPlot 
    (or it can optionally be a Spectrum type, for which a SPlot will be
    created). The splot will be rendered and the user clicks on the 
    figure. When finished, return to the terminal prompt. A dictionary is
    returned with two entries, *wave* and *data*, representing the x-y
    locations selected by the user. This can always be retrieved later by
    accessing the module member *Profile.selected*.
	
	While the user makes selections, temporary markers appear on the figure
	indicating the data point that was just selected. If a mark does not 
	appear, try moving the curser slightly and trying again. Even if the line
	goes through that point, there might not actually be data there.

<a name=FitLoc></a>
- **Fit** ( *splot*, *function* = InvertedLorentzian, *params* = None)

    Given a *splot* of type Plot.SPlot, the user selects four points on the 
    spectrum and a parameterized function is fit (an inverted Lorentzian by
    default). Optionally, *splot* can be of type spectrum and a basic SPlot
    will be created for you. If the user gives an alternative *function*, 
    *params* (parameters) must be provided. *params* is to be the first guess, 
    *p0* given to scipy...curve_fit; the user can provide them expicitely, 
    or in the form of functions with the templates `f(xarray, yarray)`
    where `xarray` and `yarray` are the `wave` and `data` arrays extracted
    between the two inner points selected by the user.
	
	*InvertedLorentian* is defined in SLiPy.Algorithms.Functions. The user does
	not need to provide *params* for the default behavior.
	
	**Example:**
	```python
	# In this example I use an alternative function simply as an illustration.
	# I show the download, calibration, and fitting procedure involved ...
	
	from slipy import Fits, Simbad, Plot, Velocity, Telluric, Profile
	from slipy import Spectrum, Observatory, SlipyError, Units as u
	from slipy.Data import Elodie
	
	# See the Data package for more information about this section
	archive = Elodie.Archive()
	
	# the one with highest S/N
	sci_file = archive.files['HD332329']
	
	# all files are in `data` member, first column is the file name 
	cal_files = [ entry[0] for entry in archive.data['HD087901A'] ]
	
	Elodie.Download(
			sci_file,                      # file name to download
			resample = (5850, 5950, 0.01), # pipeline command
			outpath = 'Science-Files'      # put under this directory
		)
		
	Elodie.Download(
			*cal_files,                    # HD087901A is Regulus
			resample = (5850, 5950, 0.01), # pipeline command
			outpath = 'Calibration-Files'  # put under this directory
		)
	
	sci_file  = Fits.Find('Science-Files')[0]  # one file
	cal_files = Fits.Find('Calibration-Files') # six files (list)
	
	spectrum = Spectrum(sci_file)       # create spectrum object
	regulus  = Fits.GetData(*cal_files) # return list of spectrum objects
	
	# remove absorption due to Earth's atmosphere
	Telluric.Correct(spectrum, *regulus)
	
	# attach `ra`, `dec` and `jd` to spectrum for HelioCorrect()
	# Simbad returns values in units of decimal degrees
	ra, dec = Simbad.Position('HD332329')
	
	# elements returned from Fits.Header() are unitless. Here we take
	# the MJD at beginning of exposure and solve for JD at middle of exposure
	MJD     = Fits.Header(sci_file, 'MJD-OBS')
	EXPTIME = Fits.Header(sci_file, 'EXPTIME')
	JD      = MJD + 2400000.5 + 0.5 * EXPTIME / 86400
	 
	 # attach to spectrum object as members
	 spectrum.ra  = ra
	 spectrum.dec = dec
	 spectrum.jd  = JD * u.day
	 
	 # Elodie as at the Observatoire de Haute-Provence, France (OHP).
	 Velocity.HelioCorrect( Observatory.OHP(), spectrum )
	 
	 # create figure
	 fig = Plot.SPlot(spectrum, label='HD332329', marker='k-')
	 fig.xlabel('Wavelength (Angstrom)', labelpad=20)
	 fig.ylabel('Normalized Flux', labelpad=20)
	 fig.xlim(5885, 5905)
	 fig.legend(frameon=False)
	 
	 # Now we need to define some parameter functions to pass to Profile.Fit()
	 # with our user function
	 from slipy.Algorithms.Functions import InvertedGaussian
	 
	 # best guess for amplitude of gaussian
	 def A(x, y):
	 	return 1 - y.min().value
		
	# best guess for mean of gaussian
	def mu(x, y):
		return x[ y.argmin() ].value
		
	# best guess for standard deviation of gaussian
	def sigma(x, y):
		return y.std().value
	
	# call the Profile.Fit() function with our user defined parameterization
	line = Profile.Fit(fig, function=InvertedGaussian, params=[A, mu, sigma])
	
	# Please select four points identifying the spectral line.
	# Outer points mark the domain of the line.
	# Inner points mark the sample of the line to fit.
	# Press <Return> after making your selections ...
	
	# now `line` is a Spectrum object generated by evaluating the 
	# `InvertedGaussian` function on the larger domain selected by the
	# user with coefficients optimized using Profile.Fit().
	# You might use it to deblend this inner line with the overlapping ones.
	
	# save the figure ...
	fig.tight_layout()
	fig.save('HD332329.png')
	```
	
    ![example](Figures/Profile-Fit.png)

    **Figure 2:** The above figure was generated by running the above code 
	snippet.
	
<br>
##<a name=MontageLoc></a>[Montage](SLiPy/Montage.py)

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

    <br>
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

    <br>
    The available filter band for each survey are as follows

    | Survey  | Bands                                                   |
    |---------|---------------------------------------------------------|
    | 'DSS'   | 'J', 'H', 'K'                                           |
    | 'SDSS'  | 'U', 'G', 'R', 'I', 'Z'                                 |
    | '2MASS' | 'DSS1B', 'DSS1R', 'DSS2B', 'DSS2R', 'DSS2IR', 'Quick-V' |

    <br>
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

    <br>
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

    <br>
    The available filter band for each survey are as follows

    | Survey  | Bands                                                   |
    |---------|---------------------------------------------------------|
    | 'DSS'   | 'J', 'H', 'K'                                           |
    | 'SDSS'  | 'U', 'G', 'R', 'I', 'Z'                                 |
    | '2MASS' | 'DSS1B', 'DSS1R', 'DSS2B', 'DSS2R', 'DSS2IR', 'Quick-V' |

    <br>
    All the member functions are the same name as in *SubField*, but now with
    the addition of a final step:

    - *Finalize* ( *resolution*, \*\**kwargs* ):

        Collect all SubField/master mosaics into a single folder and
        run Mosaic() on them for a single final image.


<br>
##<a name=ElodieLoc></a>[Elodie](Data/Elodie.py)

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

    <br>
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

    <br>
    **Example:**
    ```python
    # all files in the archive for Altair (file name is first element in pair)
    files = [ x[0] for x in archive.data['HD187642'] ]

    # download files to current directory, resample spectrum to wavelengths
    # between 5850 and 5950 with a resolution of 0.01 Angstroms per pixel
    Elodie.Download( *files, resample=(5850, 5950, 0.01) )
    ```
