# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/DataType.py 
"""
Class object representations for astronomical data.
"""
import numpy as np 
from scipy.interpolate import interp1d
from copy import deepcopy
from astropy.io import fits as pyfits
from ..Framework.Options import Options, OptionsError

class DataError(Exception):
	"""
	Exception specific to Python.Astro.DataType's.
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
	def __init__(self, argument, **kwargs ):
		"""
		Hold spectrum `data` from file name `argument`. Alternatively,
		construct spectra with a one-dimensional numpy.ndarray as `argument`.
		`wavecal` is assumed to be true for file `input`.
		"""
		try:

			self.options = Options( kwargs,
				{
					'wavecal': True     , # fit wavelength vector to data
					'crpix1' : 'crpix1' , # reference pixel header keyword 
					'crval1' : 'crval1' , # value at reference pixel 
					'cdelt1' : 'cdelt1' , # resolution (delta lambda)
				})

			self.wavecal = self.options('wavecal')
			self.crpix1  = self.options('crpix1')
			self.crval1  = self.options('crval1')
			self.cdelt1  = self.options('cdelt1')

			# observables (needed for corrections)
			self.ra  = None # right ascension
			self.dec = None # declination
			self.jd  = None # julian date

			if type(argument) is str:
				# assume we are importing from a file
				self.data = pyfits.getdata(argument)

				if self.wavecal:
					# attempt to build vector from header info
					with pyfits.open(argument) as hdulist:
						self.rpix = hdulist[0].header[self.crpix1]
						self.rval = hdulist[0].header[self.crval1]
						self.delt = hdulist[0].header[self.cdelt1]
						self.wave = WaveVector(self.rpix, self.rval,
							self.delt, np.shape(self.data)[0])

				else:
					# wave vector will just be indices of self.data
					self.wave = np.arange(len(self.data))

			elif ( type(argument) is np.ndarray and
				len(np.shape(argument)) == 1 ):
				# we are initializing from a numpy array
				self.data = argument
				# wavecal ignored on numpy.array construction
				self.wave = np.arange(len(argument))
			
			else:
				# we got the wrong argument type
				raise DataError('Spectrum() object expected either a '
				'1D numpy.ndarray or file name as a constructor argument.')
	
		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise DataError('Failed to construct Spectrum() object.')

		except IOError as err:
			print(' --> IOError:', err)
			raise DataError('Failed to construct Spectrum() object.')

	def resample(self, first, last, npix, **kwargs):
		"""
		Resample onto new wavelength pixel space. Built with numpy.linspace 
		using `first`, `last`, and `npix` as arguments.

		kwargs = {
				kind : 'linear' # same argument from scipy.interpolate.interp1d
			}
		"""
		try:
			# default function parameters
			options = Options( kwargs,
				{
					'kind':'linear' # same argument from scipy.interpolate.interp1d
				})
			# assign function parameters
			kind = options('kind')

			# check function arguments
			if first < self.wave[0] or last > self.wave[-1]:
				raise DataError('Spectrum.resample() cannot interplate outside '
				'the original domain of the flux data.')
			if npix < 1:
				raise DataError('Spectrum.resample() expects a positive '
				'integer for `npix` argument.')

			# build interplant
			f = interp1d(self.wave, self.data, kind=kind)
			# build new wave vector
			self.wave = np.linspace(first, last, npix)
			# resample data 
			self.data = f(self.wave)

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise DataError('Failed keyword assignment from '
			'Spectrum.resample().')

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
		operand.resample(self.wave[0], self.wave[-1], len(self.wave))
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
			result.data -= operand.data 
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

	def __div__(self, other):
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

	def __rdiv__(self, other):
		"""
		Right-hand division
		"""
		result = self.copy()
		result.data = other / result.data
		return result 


