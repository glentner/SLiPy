# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# SLiPy/SLiPy/DataType.py
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

class DataTypeError(SlipyError):
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
	Spectrum objects consist of a `data` vector, optionally match a
	`wavelength` vector (accessed with .data and .wave respectively).
	"""
	def __init__(self, filename, wavelengths=None, **kwargs ):
		"""
		Hold spectrum `data` from file name `argument`. Alternatively,
		construct spectra with a one-dimensional numpy.ndarray as `argument`.
		`wavecal` is assumed to be true for file `input`.
		"""
		try:

			self.options = Options( kwargs,
				{
					'wavecal': True       , # fit wavelength vector to data
					'crpix1' : 'crpix1'   , # reference pixel header keyword
					'crval1' : 'crval1'   , # value at reference pixel
					'cdelt1' : 'cdelt1'   , # resolution (delta lambda)
					'xunits' : 'Angstrom' , # units of wavelength from header
					'yunits' : 'erg cm-2 s-1'  # units of data
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

			if type(filename) is str:
				# assume we are importing from a file
				self.data = pyfits.getdata(filename)

				if self.wavecal:

					# attempt to build vector from header info
					with pyfits.open(filename) as hdulist:

						self.wave = WaveVector(
							hdulist[0].header[self.crpix1],
							hdulist[0].header[self.crval1],
							hdulist[0].header[self.cdelt1],
							np.shape(self.data)[0])

				else:

					# wave vector will just be indices of self.data
					self.wave = np.arange(len(self.data)) * u.pixel

			elif (
				len(np.shape(filename)) == 1 and
				len(np.shape(wavelengths)) == 1 and
				len(filename) == len(wavelengths)):

				# we are initializing from a numpy array
				self.data = filename.copy()
				# wavecal ignored on numpy.array construction
				self.wave = wavelengths.copy()

			else:

				# we got the wrong argument types
				raise DataTypeError('Spectrum() object expected either a set of '
				'1D numpy.ndarray`s arrays of equal length or file name as a '
				'constructor argument.')

			# add units if necessary
			if not hasattr(self.data, 'unit'):
				self.data *= u.Unit(self.yunits)

			if not hasattr(self.wave, 'unit'):
				self.wave *= u.Unit(self.xunits)

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise DataTypeError('Failed to construct Spectrum() object.')

		except IOError as err:
			print(' --> IOError:', err)
			raise DataTypeError('Failed to construct Spectrum() object.')

	def resample(self, *args, **kwargs):
		"""
		If given a single argument, it is taken to be a `Spectrum` object,
		and `self` is resampled onto the pixel space of the other spectrum.
		Otherwise, three arguments are expected. The first and second argument
		should define the lower and upper wavelength value of a domain,
		respectively. The third argument should be the number of elements
		(pixels) for the new domain. Think numpy.linspace().

		kwargs = {
				kind : 'linear' # passed to scipy.interpolate.interp1d
			}
		"""

		try:
			options = Options( kwargs, {

				'kind' : 'linear' # passed to scipy.interpolate.interp1d
			})

			kind = options('kind')

		except OptionsError as err:
			print(' --> OptionsError: ', err)
			raise DataTypeError('From Spectrum.resample(), failed keyword '
			'assignment!')

		# build interpolation function
		f = interp1d(self.wave, self.data, kind=kind)

		if len(args) == 1:

			spectrum = args[0]

			if type(spectrum) is not Spectrum:
				raise DataTypeError('When `resample()` is given a single '
				'argument, it is expected to by of type `Spectrum`.')

			# check domain limits
			if spectrum.wave[0] < self.wave[0] or spectrum.wave[-1] > self.wave[-1]:
				raise DataTypeError('Spectrum.resample() cannot interplate outside '
				'the original domain of the flux data.')

			# build new wave-array from spectrum
			self.wave = spectrum.wave

			# bring back to original units for f()
			self.wave = self.wave.to(u.Unit(self.xunits))

			# resample data
			self.data = f(self.wave) * u.Unit(self.yunits)

			# conver to units of the given spectrum
			self.wave = self.wave.to(u.Unit(spectrum.wave.unit))
			self.data = self.data.to(u.Unit(spectrum.data.unit))

		elif len(args) == 3:

			low, high, npix = args

			# check function arguments
			if low < self.wave[0] or high > self.wave[-1]:
				raise DataTypeError('Spectrum.resample() cannot interplate outside '
				'the original domain of the flux data.')

			if npix < 1:
				raise DataTypeError('Spectrum.resample() expects a positive '
				'integer for `npix` argument.')

			# build new wave vector
			self.wave = np.linspace(low, high, npix)

			# bring back to original units for f()
			self.wave = self.wave.to(u.Unit(self.xunits))

			# resample data
			self.data = f(self.wave) * u.Unit(self.yunits)

		else:
			raise DataTypeError('Unacceptable number of arguments passed to '
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
		operand = other.copy()
		operand.resample(self)
		return self.copy(), operand

	def __add__(self, other):
		"""
		Addition
		"""
		if type(other) is Spectrum:
			# resample the data and return `new` spectrum
			result, operand = self.__new_pair(other)
			result.data += operand.data
			return result

		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data += other
		return result

	def __sub__(self, other):
		"""
		Subtraction
		"""
		if type(other) is Spectrum:
			# resample the data and return `new` spectrum
			result, operand = self.__new_pair(other)
			result.data = result.data - operand.data
			return result

		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data -= other
		return result

	def __mul__(self, other):
		"""
		Multiplication
		"""
		if type(other) is Spectrum:
			# resample the data and return `new` spectrum
			result, operand = self.__new_pair(other)
			result.data *= operand.data
			return result

		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data *= other
		return result

	def __truediv__(self, other):
		"""
		Division
		"""
		if type(other) is Spectrum:
			# resample the data and return `new` spectrum
			result, operand = self.__new_pair(other)
			result.data /= operand.data
			return result

		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data /= other
		return result

	def __iadd__(self, other):
		"""
		In-place addition
		"""
		if type(other) is Spectrum:
			# resample operand before operation
			operand = other.copy()
			operand.resample(self.wave[0], self.wave[-1], len(self.wave))
			self.data += operand.data

		# otherwise, assume we have a scalar (be careful)
		else: self.data += other

	def __isub__(self, other):
		"""
		In-place subtraction
		"""
		if type(other) is Spectrum:
			# resample operand before operation
			operand = other.copy()
			operand.resample(self.wave[0], self.wave[-1], len(self.wave))
			self.data -= operand.data

		# otherwise, assume we have a scalar (be careful)
		else: self.data -= other

	def __imul__(self, other):
		"""
		In-place multiplication
		"""
		if type(other) is Spectrum:
			# resample operand before operation
			operand = other.copy()
			operand.resample(self.wave[0], self.wave[-1], len(self.wave))
			self.data *= operand.data

		# otherwise, assume we have a scalar (be careful)
		else: self.data *= other

	def __itruediv__(self, other):
		"""
		In-place division
		"""
		if type(other) is Spectrum:
			# resample operand before operation
			operand = other.copy()
			operand.resample(self.wave[0], self.wave[-1], len(self.wave))
			self.data /= operand.data

		# otherwise, assume we have a scalar (be careful)
		else: self.data /= other

	def __radd__(self, other):
		"""
		Right-hand Addition
		"""
		result = self.copy()
		result.data += other
		return result

	def __rsub__(self, other):
		"""
		Right-hand subtraction
		"""
		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data = other - result.data
		return result

	def __rmul__(self, other):
		"""
		Right-hand multiplication
		"""
		# otherwise, assume we have a scalar (be careful)
		result = self.copy()
		result.data *= other
		return result

	def __rtruediv__(self, other):
		"""
		Right-hand division
		"""
		result = self.copy()
		result.data = other / result.data
		return result
