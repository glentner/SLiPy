# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# slipy/SLiPy/Correlate.py
"""
Correlate - Module of correlation functions for astronomical data.
"""
import numpy as np

from .DataType import Spectrum, DataError
from ..Framework.Options import Options, OptionsError

class CorrelateError(Exception):
	"""
	Exception specific to Correlate module.
	"""
	pass

def RMS( array ):
	"""
	RMS( array ):

	Return the root mean square of an `array`.
	"""
	if type(array) is not np.ndarray or len(np.shape(array)) != 1:
		raise CorrelateError('RMS() expects 1D numpy.ndarray`s.')

	return np.sqrt( ( array**2 ).sum() ) / len(array)

def Xcorr( spectrumA, spectrumB, **kwargs ):
	"""
	Xcorr( spectrumA, spectrumB, **kwargs ):

	Cross correlate two spectra of equal pixel length. The function returns
	an integer value representing the best shift within a `lag` based on
	the computed RMS of each configuration.
	"""
	try:
		# define `options` for Xcorr()
		options = Options( kwargs,
			{
				'lag' : 25 # pixels to shift for correlation
			})

		# check argument types, values
		if ( type(spectrumA) is not Spectrum or
			type(spectrumB) is not Spectrum ):
			raise CorrelateError('Xcorr() expects `Spectrum` arguments.')

		elif len(spectrumA.data) != len(spectrumB.data):
			raise CorrelateError('Xcorr() expects `Spectrum` arguments to '
			'be of equal length.')

		# assign `lag` and the pixel length of the spectra
		lag  = options('lag')
		npix = len(spectrumA.data)

		# shift spectra `left` over each other
		left = np.array([ RMS(diff) for diff in [ spectrumA.data[-shift:] -
			spectrumB.data[:shift] for shift in range(-lag,0) ] ])

		# shift spectra `right` over each other
		right = np.array([ RMS(diff) for diff in [ spectrumA.data[:-shift] -
			spectrumB.data[shift:] for shift in range(1,lag+1) ] ])

		# include `zero` shift in rms vector.
		rms = np.hstack((left, RMS(spectrumA.data - spectrumB.data), right))

		# return the shift corresponding to the minimum RMS
		return rms.argmin() - lag

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise CorrelateError('Inappropiate keyword arguments passed '
		' to Xcorr().')
