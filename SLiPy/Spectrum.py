# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# SLiPy/SLiPy/Spectrum.py
"""
Class object representations for astronomical data.
"""
import numpy as np
from scipy.interpolate import interp1d
from copy import deepcopy
from astropy.io import fits as pyfits
from astropy import units as u

from .. import SlipyError
from ..Framework.Options import Options, OptionsError

class SpectrumError(SlipyError):
	"""
	Exception specific to this module.
	"""
	pass

def WaveVector( rpix, rval, delt, npix ):
	"""
	WaveVector(rpix, rval, delt, npix):

	Construct numpy array of wavelength values based on:

		`rpix` : reference pixel index
		`rval` : wavelength at reference pixel
		`delt` : resolutions (delta lambda)
		`npix` : length of desired array
	"""
	if rpix < 0:
		raise FitsError('`rpix` must be a positive integer '
		'for WaveVector().'.format(rpix))

	if rval <= 0:
		raise FitsError('`rval` wavelength should be a positive '
		'value for WaveVector().'.format(rval))

	if delt <= 0:
		raise FitsError('`delt` wavelength should be a positive '
		'value for Wavevector().'.format(delt))

	if npix < 1:
		raise FitsError('`npix` must be a positive integer '
		'for WaveVector().'.format(npix))

	return rval + ( np.array(range(npix)) - rpix + 1 ) * delt

class Spectrum:
    """
    The Spectrum object is a container class for a data array and its
    corresponding wavelength array. See the __init__ documentation for
    more information on constructing spectrum objects. 
    
    Addition, subtraction, multiplication, and division (including in-place 
    operations, e.g., '+=') are defined for both spectrum on spectrum 
    operations as well as scalar operations. When two spectra are operated on, 
    the LHS spectrum takes precedent. One of the spectra must be contained 
    within the other (i.e., their wavelength domain is either equal to or 
    entirely interior to the other). The outer spectrum is resampled onto the 
    inner spectrum's pixel space and the operation is applied element-wise. 
    To state it briefly, only the affected region of the LHS spectrum is 
    operated on. This makes units dangerous for multiplication and division.
    Only in the case of multiplying/dividing spectra of equivilent wavelength
    arrays will the units be properly applied. Otherwise, the RHS units will
    be ignored. Considering the physical context within which it makes sense
    to apply these operations to spectra this is justified; the data will have
    `dimensionless` units in all likelihood. For scalar operations, units are
    implied to be the same as the data for addition/subtraction.
    
    The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
    to return a binary Spectrum of True and False values. The same convention
    applies as above. Either the LHS or RHS must be contained within the other,
    and the LHS is compared on the overlapping regions. Data outside this 
    overlapping region is returned as False.
    
    The shift operations (>> and <<) are defined to mean a shift in the
    spectrum. The addition and subtraction operate on the `data`. These
    operations apply to the `wave` array. `Spectrum << 2 * u.Angstrom` say
    to blue shift the spectrum by 2 Angstroms. If the RHS is a dimensionless
    number, the wavelength units are implied. Also, one can shift a spectrum
    by another spectrum (e.g., `spectrumA >> spectrumB`), where the `wave` 
    arrays would be operated on only. In these cases, they should
    have the same number of pixels! Finally, to get a variable shift across
    the spectrum without creating a whole spectrum with junk data, you can
    shift using a numpy array (also of equal length). I no units are detected,
    the wavelength units will be implied.
    
    Other operations:
    
        SpectrumA in SpectrumB 
        
            --> The wavelength domain of A is either equal to or contained 
                by B. True/False
        
        len(Spectrum)
        
            --> Number of pixels.
        
        str(Spectrum)
        
            --> Calls str() on member arrays
        
        []
        
            --> Access the data via indexing. Either access a single
                pixel or a range (e.g., Spectrum[3], Spectrum[4:-5]).
                Alternatively, given a Quantity, return an estimate of
                the value of the spectrum at that location (e.g., 
                Spectrum[588.89 * u.nm]) will return a linear 
                approximation at that wavelength location.
        
    Member Functions:
    
        .resample(*args, kind='linear')
        
            Resample spectrum given new domain. Arguments can either be
            a different spectrum, or a set of three arguments passed to
            numpy.linspace().
        
        .insert(spectrum, kind='linear')
        
            All pixels in `self` interior to the `spectrum` are replaced
            with a resampling based on the `spectrum`.
        
        .copy()
        
            Calls copy.deepcopy() on `self`. To say SpectrumA = SpectrumB
            is to say the SpectrumA *is* SpectrumB. Replace with 
            SpectrumA = SpectrumB.copy()
            
    """
    def __init__(self, *args, **kwargs ):
        """
        There are a few ways to create a Spectrum. If a single string 
        argument is given, it is taken as a file name and used to read in
        data from a FITS file. With the keyword argument `wavecal` set as
        True (the default case), elements are read from the header to create
        a corresponding wavelength array to go with the data.
        
        If an array-like object is given, it is converted to a numpy array and
        taken as the spectrum data. A wavelength array will be generated that 
        is simply an index (pixel) count. But if a second argument is provided, 
        it will serve as the wavelength array. These must be equal in length 
        however.
        
        Units will be imposed for these arrays. When initialized from a file,
        the default units are `Angstrom` and `dimensionless_unscaled` for the
        wavelength and data arrays respectively. Alternatives can be applied
        by providing the keyword arguments `xunits` and `yunits`. If 
        initialized via an array-like object, `dimensionless` will
        only be applied if no units are detected. If no wavelength array is
        provided, the generated wavelength array will have `pixel` units.
        Units are again only applied if none are detected for the given array.
        """
        try:

            self.options = Options( kwargs,
                {
                    'wavecal': True       , # fit wavelength vector to data
                    'crpix1' : 'crpix1'   , # reference pixel header keyword
                    'crval1' : 'crval1'   , # value at reference pixel
                    'cdelt1' : 'cdelt1'   , # resolution (delta lambda)
                    'xunits' : 'Angstrom' , # units of wavelength from header
                    'yunits' : ''           # units of data
                })

            self.wavecal = self.options('wavecal')
            self.crpix1  = self.options('crpix1')
            self.crval1  = self.options('crval1')
            self.cdelt1  = self.options('cdelt1')
            self.xunits  = self.options('xunits')
            self.yunits  = self.options('yunits')

            # observables (needed for corrections)
            self.ra  = None # right ascension
            self.dec = None # declination
            self.jd  = None # julian date

            if len(args) == 1:
            
                if type(args[0]) is str:
                
                    # assume we are importing from a file
                    filename = args[0]
                    
                    self.data = pyfits.getdata(filename)

                    if self.wavecal:

                        # attempt to build vector from header info
                        with pyfits.open(filename) as hdulist:

                            self.wave = WaveVector(
                                hdulist[0].header[self.crpix1],
                                hdulist[0].header[self.crval1],
                                hdulist[0].header[self.cdelt1],
                                np.shape(self.data)[0]
                                ) * u.Unit(self.xunits)

                    else:

                        # wavelength array will just be pixel numbers 
                        self.wave = ( np.arange(len(self.data), dtype = np.int_) 
                            + 1) * u.pixel
                    
                elif hasattr(args[0], '__iter__'):
                    
                    # we have an array-like object, use forgiving procedure
                    # to create numpy array
                    data = args[0]
                    
                    if hasattr(data[0], 'unit'):
                        self.data = np.array([x.value for x in data], 
                            dtype = np.float_) * data[0].unit
                    
                    else:
                        self.data = np.array(list(data), 
                            dtype = np.float_) * u.Unit(self.yunits)
                    
                    # wavelength array will just be pixel numbers
                    self.wave = ( np.arange(len(data), dtype = np.int_ ) 
                        + 1) * u.pixel
                    
                    # create copy of data to disconnect references
                    self.data = self.data.copy()
                    
                else:
                    raise SpectrumError('Spectrum objects can only be '
                    'initialized by a file name or an array-like object!')

            elif len(args) == 2:

                # arguments should be two array-like objects
                data, wave = args
                
                if not hasattr(data,'__iter__') or not hasattr(wave,'__iter__'):
                    raise SpectrumError('Spectrum object expects arguments to '
                    'be array-like when two are given.')
                
                if len(data) != len(wave):
                    raise SpectrumError('Spectrum object expects array '
                    'array arguments to be of equal length.')
                    
                if hasattr(data[0], 'unit'):
                    self.data = np.array([x.value for x in data], 
                        dtype = np.float_) * data[0].unit
                
                else:
                    self.data = np.array(list(data), 
                        dtype = np.float_) * u.Unit(self.yunits)
                        
                if hasattr(wave[0], 'unit'):
                    self.wave = np.array([x.value for x in wave], 
                        dtype = np.float_) * wave[0].unit
                
                else:
                    self.wave = np.array(list(wave), 
                        dtype = np.float_) * u.Unit(self.xunits)
                        
                # create copy of data to disconnect references
                self.data = self.data.copy()
                self.wave = self.wave.copy()

            else:

                # we got the wrong argument number
                raise SpectrumError('Spectrum objects must be initialized '
                'with either one or two arguments.')

        except OptionsError as err:
            print(' --> OptionsError:', err)
            raise SpectrumError('Failed to construct Spectrum() object.')

        except IOError as err:
            print(' --> IOError:', err)
            raise SpectrumError('Failed to construct Spectrum() object.')

    def resample(self, *args, kind = 'linear'):
        """
        If given a single argument, it is taken to be a `Spectrum` object,
        and `self` is resampled onto the pixel space of the other spectrum.
        Otherwise, three arguments are expected. The first and second argument
        should define the lower and upper wavelength value of a domain,
        respectively. The third argument should be the number of elements
        (pixels) for the new domain. Think numpy.linspace().
        """

        # build interpolation function
        f = interp1d(self.wave, self.data, kind=kind)

        if len(args) == 1:

            spectrum = args[0]

            if type(spectrum) is not Spectrum:
                raise SpectrumError('When `resample()` is given a single '
                'argument, it is expected to by of type `Spectrum`.')

            # check domain limits
            if spectrum not in self:
                raise SpectrumError('Spectrum.resample() cannot interplate outside '
                'the original domain of the flux data.')

            # build new wave-array from spectrum
            self.wave = spectrum.wave.to(self.wave.unit)

            # resample data
            self.data = f(self.wave) * self.data.unit

            # convert to units of the given spectrum
            self.wave = self.wave.to(u.Unit(spectrum.wave.unit))
            self.data = self.data.to(u.Unit(spectrum.data.unit))

        elif len(args) == 3:

            low, high, npix = args
            
            if not hasattr(low, 'unit'):
                low *= self.wave.unit
            
            if not hasattr(high, 'unit'):
                high *= self.wave.unit
            
            if hasattr(npix, 'unit'):
                if npix.unit != u.pixel:
                    raise SpectrumError('`npix` can only have units of `pixels` '
                    'if given units!')
                npix = npix.value

            # check function arguments
            if low < self.wave[0] or high > self.wave[-1]:
                raise SpectrumError('Spectrum.resample() cannot interplate outside '
                'the original domain of the flux data.')

            if npix < 2:
                raise SpectrumError('Spectrum.resample() expects an integer '
                'greater than 2 for `npix` argument.')

            # build new wave vector
            if not hasattr(low, 'unit'): low *= self.wave.unit
            if not hasattr(high, 'unit'): high *= self.wave.unit
            self.wave = np.linspace(low, high, npix).to(self.wave.unit)

            # resample data
            self.data = f(self.wave) * self.data.unit

        else:
            raise SpectrumError('Unacceptable number of arguments passed to '
            'Spectrum.resample()')

    def copy(self):
        """
        Call to copy.deepcopy on `self`.
        """
        return deepcopy(self)

    def __new_pair(self, other):
        """
        Create `copy`s of `self` and `other` for operator overloading.
        """
        operand, result = other.copy(), self.copy()

        if operand in result:
            result.resample(operand)

        elif result in operand:
            operand.resample(result)

        else:
            raise SpectrumError('These spectra do not overlap!')

        return result, operand

    def __contains__(self, other):
        """
        Check whether the domain of `other` is contained within the
        domain of `self`.
        """
        if other.wave[0] >= self.wave[0] and other.wave[-1] <= self.wave[-1]:
            return True
        
        else:
            return False
    
    def insert(self, other, kind = 'linear'):
        """
        Given a Spectrum, `other`, contained within the wavelength domain
        of `self`, replace all pixels in the overlapping region with that
        of an interpolation built on `other`. `kind` is passed interp1d.
        """
        if type(other) is not Spectrum:
            raise SpectrumError('Spectrum.insert() expects an argument of '
            'type Spectrum!')
        
        if other not in self:
            raise SpectrumError('The domain of the `other` spectrum is not '
            'contained within this spectrum!')
        
        # interpolant of the `other` spectrum
        f = interp1d(other.wave.to(self.wave.unit), other.data, kind = kind)
        
        self.data = np.hstack((
            
            # the lhs (not affected)
            self.data[np.where(self.wave < other.wave[0])],
            
            # the interpolated pixels covered by `other`
            f(self.wave[np.where(np.logical_and(self.wave >= other.wave[0], 
                self.wave <= other.wave[-1]))]) * self.data.unit,
                
            # the rhs (not affected)
            self.data[np.where(self.wave > other.wave[-1])]
        
        # the units get messed up, so we have to re-initialize them
        )).value * self.data.unit

    def __str__(self):
        """
        Printing a Spectrum displays the two arrays.
        """
        return str(self.data) + '\n' + str(self.wave)
    
    def __repr__(self):
        """
        Same as __str__.
        """
        return str(self)
    
    def __len__(self):
        """
        Number of pixels (resolution) in spectrum.
        """
        return len(self.wave)
        
    def __add__(self, other):
        """
        Addition. When two spectrum objects are added, one needs
        to be entirely contained by the other (or at least have the same
        domain). The overlapping regions are extracted, operated on element
        wise, and re-inserted into the lhs spectrum. For scalar operations, 
        the operation is performed element wise. If no units are given for
        the scalar addition, it is assumed that the units are the same as
        for the data. Be careful, if two Spectrum objects have incompatible
        units for their data, this will raise an error.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum addition can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data + rhs.data).to(self.data.unit)
            
            result = self.copy()
            result.insert(lhs)
            
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = result.data + other
        return result

    def __sub__(self, other):
        """
        Subtraction. When two spectrum objects are subtracted, one needs
        to be entirely contained by the other (or at least have the same
        domain). The overlapping regions are extracted, operated on element
        wise, and re-inserted into the lhs spectrum. For scalar operations, 
        the operation is performed element wise. If no units are given for
        the scalar subtraction, it is assumed that the units are the same as
        for the data. Be careful, if two Spectrum objects have incompatible
        units for their data, this will raise an error.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum subtraction can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data - rhs.data).to(self.data.unit)
            
            result = self.copy()
            result.insert(lhs)
            
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = result.data - other
        return result

    def __mul__(self, other):
        """
        Multiplication. When two spectrum objects are multiplied, one needs
        to be entirely contained by the other (or at least have the same
        domain). The overlapping regions are extracted, operated on element
        wise, and re-inserted into the lhs spectrum. For scalar operations, 
        the operation is performed element wise. Be careful, if two Spectrum 
        objects have incompatible units for their data, this will raise an error.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum mulitplication can only work '
                'with pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = lhs.data * rhs.data
            
            result = self.copy()
            result.insert(lhs)
            
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        result.data = result.data * other
        return result

    def __truediv__(self, other):
        """
        Division. When two spectrum objects are divided, one needs
        to be entirely contained by the other (or at least have the same
        domain). The overlapping regions are extracted, operated on element
        wise, and re-inserted into the lhs spectrum. For scalar operations, 
        the operation is performed element wise. Be careful, if two Spectrum 
        objects have incompatible units for their data, this will raise an error.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum division can only work '
                'with pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = lhs.data / rhs.data
            
            result = self.copy()
            result.insert(lhs)
            
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        result.data = result.data / other
        return result

    def __iadd__(self, other):
        """
        In-place addition.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum addition can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data + rhs.data).to(self.data.unit)
            
            self.insert(lhs)
            
            return self

        # otherwise, assume we have a scalar (be careful)
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= self.data.unit
        self.data = self.data + other
        return self

    def __isub__(self, other):
        """
        In-place subtraction.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum subtraction can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data + rhs.data).to(self.data.unit)
            
            self.insert(lhs)
            
            return self

        # otherwise, assume we have a scalar (be careful)
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= self.data.unit
        self.data = self.data - other
        return self

    def __imul__(self, other):
        """
        In-place multiplication.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum multiplication can only work '
                'with pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = lhs.data * rhs.data
            
            self.insert(lhs)
            
            return self

        # otherwise, assume we have a scalar (be careful)
        # if no units are detected, implicitely assume data units
        self.data = self.data * other
        return self

    def __itruediv__(self, other):
        """
        In-place division.
        """
        
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum division can only work '
                'with pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = lhs.data / rhs.data
            
            self.insert(lhs)
            
            return self

        # otherwise, assume we have a scalar (be careful)
        # if no units are detected, implicitely assume data units
        self.data = self.data / other
        return self

    def __radd__(self, other):
        """
        Right-hand addition.
        """
        # this only makes sense with scalar operations
        result = self.copy()
        if not hasattr(other, 'unit'): other *= self.data.unit
        result.data = other + result.data
        return result

    def __rsub__(self, other):
        """
        Right-hand subtraction.
        """
        # this only makes sense with scalar operations
        result = self.copy()
        if not hasattr(other, 'unit'): other *= self.data.unit
        result.data = other - result.data
        return result

    def __rmul__(self, other):
        """
        Right-hand multiplication.
        """
        # this only makes sense with scalar operations
        result = self.copy()
        result.data = other * result.data
        return result

    def __rtruediv__(self, other):
        """
        Right-hand division.
        """
        # this only makes sense with scalar operations
        result = self.copy()
        result.data = other / result.data
        return result
    
    def __lshift__(self, other):
        """
        Left shift operator.
        
        The shift operations (>> and <<) are defined to mean a shift in the
        spectrum. The addition and subtraction operate on the `data`. These
        operations apply to the `wave` array. `Spectrum << 2 * u.Angstrom` say
        to blue shift the spectrum by 2 Angstroms. If the RHS is a dimensionless
        number, the wavelength units are implied. Also, one can shift a spectrum
        by another spectrum (e.g., `spectrumA >> spectrumB`), where the `wave` 
        arrays would be operated on only. In these cases, they should
        have the same number of pixels! Finally, to get a variable shift across
        the spectrum without creating a whole spectrum with junk data, you can
        shift using a numpy array (also of equal length). I no units are detected,
        the wavelength units will be implied.
        """
        
        if type(other) is Spectrum:
            
            if len(self) != len(other):
                raise SpectrumError('The left shift operator requires the '
                'length of both spectrum objects be the same!')
            
            result = self.copy()
            result.wave = result.wave - other.wave
            return result
        
        elif hasattr(other, 'unit'):
            
            if hasattr(other, '__init__'):
                if len(self) != len(other):
                    raise SpectrumError('The left shift operator requires '
                    'the length of both the spectrum and the RHS to have the '
                    'same length!')
            
            result = self.copy()
            result.wave = result.wave - other
            return result
        
        elif type(other) is np.ndarray:
            
            if len(self) != len(other):
                raise SpectrumError('The left shift operator requires the '
                'length of both the spectrum and the RHS to have the '
                'same length!')
            
            result = self.copy()
            result.wave = result.wave - other * result.wave.unit
            return result
        
        else:
            
            if not hasattr(other, 'unit'):
                other *= self.wave.unit
            
            result = self.copy()
            result.wave = result.wave - other
            return result
            
    def __rshift__(self, other):
        """
        Right shift operator.
        
        The shift operations (>> and <<) are defined to mean a shift in the
        spectrum. The addition and subtraction operate on the `data`. These
        operations apply to the `wave` array. `Spectrum << 2 * u.Angstrom` say
        to blue shift the spectrum by 2 Angstroms. If the RHS is a dimensionless
        number, the wavelength units are implied. Also, one can shift a spectrum
        by another spectrum (e.g., `spectrumA >> spectrumB`), where the `wave` 
        arrays would be operated on only. In these cases, they should
        have the same number of pixels! Finally, to get a variable shift across
        the spectrum without creating a whole spectrum with junk data, you can
        shift using a numpy array (also of equal length). I no units are detected,
        the wavelength units will be implied.
        """
        
        if type(other) is Spectrum:
            
            if len(self) != len(other):
                raise SpectrumError('The right shift operator requires the '
                'length of both spectrum objects be the same!')
            
            result = self.copy()
            result.wave = result.wave + other.wave
            return result
        
        elif hasattr(other, 'unit'):
            
            if hasattr(other, '__init__'):
                if len(self) != len(other):
                    raise SpectrumError('The right shift operator requires '
                    'the length of both the spectrum and the RHS to have the '
                    'same length!')
            
            result = self.copy()
            result.wave = result.wave + other
            return result
        
        elif type(other) is np.ndarray:
            
            if len(self) != len(other):
                raise SpectrumError('The right shift operator requires the '
                'length of both the spectrum and the RHS to have the '
                'same length!')
            
            result = self.copy()
            result.wave = result.wave + other * result.wave.unit
            return result
        
        else:
            
            if not hasattr(other, 'unit'):
                other *= self.wave.unit
            
            result = self.copy()
            result.wave = result.wave + other
            return result

    def __lt__(self, other):
        """
        Less than, '<'
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data < rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data < other) * u.dimensionless_unscaled
        return result
        
    def __le__(self, other):
        """
        Less than or equal to, '<='
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data <= rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data <= other) * u.dimensionless_unscaled
        return result
        
    def __gt__(self, other):
        """
        Greater than, '>'
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data > rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data > other) * u.dimensionless_unscaled
        return result

    def __ge__(self, other):
        """
        Greater than or equal to, '>='
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data >= rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data >= other) * u.dimensionless_unscaled
        return result
        
    def __eq__(self, other):
        """
        Equal to, '=='
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data == rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data == other) * u.dimensionless_unscaled
        return result

    def __ne__(self, other):
        """
        Not equal to, '!='
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = (lhs.data != rhs.data) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = (result.data != other) * u.dimensionless_unscaled
        return result
        
    def __and__(self, other):
        """
        Logical and, '&'
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = np.logical_and(lhs.data.value, 
                rhs.data.to(lhs.data.unit).value) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_and(result.data.value, 
            other.to(result.data.unit).value) * u.dimensionless_unscaled
        return result

    def __or__(self, other):
        """
        Logical or, '|'
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = np.logical_or(lhs.data.value, 
                rhs.data.to(lhs.data.unit).value) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_or(result.data.value, 
            other.to(result.data.unit).value) * u.dimensionless_unscaled
        return result

    def __xor__(self, other):
        """
        Logical, exclusive or, '^'
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """
        if type(other) is Spectrum:
            
            if ((self.wave.unit == u.pixel or other.wave.unit == u.pixel) and
                len(self) != len(other) ):
                raise SpectrumError('Spectrum comparisons can only work with '
                'pixel units when both wavelength arrays have pixel '
                'units and are of the same length!')
            
            # resample the data and return `new` spectrum objects
            lhs, rhs = self.__new_pair(other)
            lhs.data = np.logical_xor(lhs.data.value, 
                rhs.data.to(lhs.data.unit).value) * u.dimensionless_unscaled
            
            result = self.copy()
            
            # assume false
            result.data = np.zeros(len(self)) * u.dimensionless_unscaled
            result.insert(lhs)
            return result

        # otherwise, assume we have a scalar (be careful)
        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_xor(result.data.value, 
            other.to(result.data.unit).value) * u.dimensionless_unscaled
        return result
        
    def __rand__(self, other):
        """
        Logical and, '&' - Right Handed (scalar)
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """

        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_and(other.to(result.data.unit).value, 
            result.data.value) * u.dimensionless_unscaled
        
        return result

    def __ror__(self, other):
        """
        Logical or, '|' - Right Handed (scalar)
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """

        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_or(other.to(result.data.unit).value, 
            result.data.value) * u.dimensionless_unscaled
        
        return result

    def __rxor__(self, other):
        """
        Logical, exclusive or, '^' - Right Handed (scalar)
        
        The comparison operations ( >, <, >=, <=, ==, !=, &, ^, | ) are defined
        to return a binary Spectrum of True and False values. The same convention
        applies as above. Either the LHS or RHS must be contained within the other,
        and the LHS is compared on the overlapping regions. Data outside this 
        overlapping region is returned as False.
        """

        result = self.copy()
        
        # if no units are detected, implicitely assume data units
        if not hasattr(other, 'unit'): other *= result.data.unit
        
        result.data = np.logical_xor(other.to(result.data.unit).value, 
            result.data.value) * u.dimensionless_unscaled
        
        return result
    
    def __getitem__(self, key):
        """
        Indexing and Slicing. 
        
        If the `key` is a slice object, extract
        the segment of the spectrum that falls within that wavelength
        regime (units assumed to be that of the `wave` array!)
        
        Otherwise, return a linear approximation of the spectrum at that
        wavelength. If the spectrum has pixel units, and given an integer,
        the exact value is returned. If a decimal is given, return an 
        approximation between the two pixels. If the wavelength array is
        in units of length, assume the units of the array if none are given.
        """
        
        if isinstance(key, slice):
            
            start, stop, step = key.start, key.stop, key.step
            
            if not start:
                start = self.wave[0].value
            
            if not stop:
                stop = self.wave[-1].value
            
            if start < self.wave[0].value:
                raise SpectrumError('{} is outside the domain of the '
                'spectrum: {} -> {}'.format(start * self.wave.unit,
                self.wave[0], self.wave[-1]))
            
            if stop > self.wave[-1].value:
                raise SpectrumError('{} is outside the domain of the '
                'spectrum: {} -> {}'.format(stop * self.wave.unit,
                self.wave[0], self.wave[-1]))
            
            if step and (step > self.wave[-1].value - self.wave[0].value):
                raise SpectrumError('The specified `step` in the slice '
                '({}) is greater than the domain of the spectrum: '
                '{} -> {}'.format(step * self.wave.unit, self.wave[0],
                self.wave[-1]))
            
            if not step:
                
                # we are just going to return a segment of the spectrum
                # on its own pixel space
                
                return Spectrum(
                    
                    self.data.value[np.where(np.logical_and(
                            
                            self.wave >= start * self.wave.unit,
                            self.wave <= stop  * self.wave.unit
                        
                        ))] * self.data.unit,
                    
                    self.wave.value[np.where(np.logical_and(
                        
                            self.wave >= start * self.wave.unit,
                            self.wave <= stop  * self.wave.unit
                        
                        ))] * self.wave.unit
                    )
            
            else:
                
                # we are going to resample the spectrum on the following
                # wavelength domain
                
                result = self.copy()
                result.resample(start, stop, int(1 + (stop - start) / step))
                
                return result
            
        else:
            
            # return the exact value, or a linear approximation of the 
            # spectrum at the location of `key`
            
            if hasattr(key, 'unit'):
                key = key.to(self.wave.unit).value
            
            if key < self.wave[0].value:
                raise SpectrumError('{} is outside the domain of the '
                'spectrum: {} -> {}'.format(key * self.wave.unit,
                self.wave[0], self.wave[-1]))
            
            if key > self.wave[-1].value:
                raise SpectrumError('{} is outside the domain of the '
                'spectrum: {} -> {}'.format(key * self.wave.unit,
                self.wave[0], self.wave[-1]))
            
            return interp1d(self.wave, self.data)(key) * self.data.unit