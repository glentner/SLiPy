# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/Astro/DataType.py 
"""
Class object representations for astronomical data.
"""
import pyfits, numpy as np 

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
	Spectra consist of a `data` vector, optionally match a  
	`wavelength` vector (accessed with .data and .wave respectively).
	"""
	def __init__(self, filename, **kwargs ):
		"""
		Hold spectrum `data` from `filename`, optionally build `wave` 
		vector.
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

			self.data = pyfits.getdata(filename)

			if self.wavecal:

				with pyfits.open(filename) as hdulist:
					self.rpix = hdulist[0].header[self.crpix1]
					self.rval = hdulist[0].header[self.crval1]
					self.delt = hdulist[0].header[self.cdelt1]
					self.wave = WaveVector( 
						self.rpix, self.rval, self.delt, np.shape(self.data)[0]
						)
				
		except OptionsError as err:
			print(' --> OptionsError:', err.msg)
			raise DataError('Failed to construct spectrum.')

		except IOError as err:
			print(' --> IOError:', err)
			raise DataError('Failed to construct spectrum.')

	def __resample(self, spectra):
		"""
		Verify wavelength vector of `spectra` for operator
		overloading. Resample `spectra` onto `self.wave` before 
		applying operator to `self.data`.
		"""
		# check that we have a WaveVector 
		if type(spectra) is not Spectrum:
			raise TypeError('Spectrum.__resample needs type Spectrum')
		elif not spectra.wavecal:
			raise DataError('Spectrum.__resample received a Spectrum'
			' without a WaveVector.')

		# quick check 
        #if ( len(spectra.wave) == len(self.wave) and

