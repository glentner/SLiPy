# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/SLiPy/Tellucic.Py
"""
Telluric - Corrections for atmospheric absoption lines.
"""
import numpy as np
from astropy import units as u

from .. import SlipyError
from ..Framework.Options import Options, OptionsError
from .Correlate import Xcorr, CorrelateError
from .DataType import Spectrum, DataTypeError

class TelluricError(SlipyError):
	"""
	Exception specific to the Telluric module.
	"""
	pass

def Correct(spectrum, *calibration, **kwargs):
	"""
	Correct(spectrum, *calibration, **kwargs):

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
	try:
		# default keyword arguments
		options = Options( kwargs,
			{
				'lag'  : 25            , # pixel shift limit for XCorr()
				'range':(0.5, 2.0, 151)  # numpy.linspace for amplitude trials
			})

		# check arguments
		if not calibration:
			raise TelluricError('At least one `calibration` spectrum '
			'needs to be provided for Telluric.Correction().')

		if type(spectrum) is not Spectrum:
			raise TelluricError('Telluric.Correct() expects all arguments '
			'of type Spectrum.')

		for cal in calibration:
			if type(cal) is not Spectrum:
				raise TelluricError('Telluric.Correct() expects all '
				'arguments to be of type Spectrum.')
			if len(spectrum.data) != len(cal.data):
				raise TelluricError('Telluric.Correct() expects all '
				'Spectrum arguments to have an equal number of elements.')

		if len(options('range')) != 3:
			raise OptionsError('`range` expects a tuple of length 3.')

		# assign parameters
		lag    = options('lag')
		amp    = np.linspace( *options('range') )
		trials = len(amp)
		npix   = len(spectrum.data)

		# quit if too big of a task
		if trials*npix > 1e8:
			raise TelluricError('Telluric.Correct() is programmed to quit '
			'if it detects a request to operate on matrices with more '
			'than 10**8 elements.')

		# resample all the calibration spectra
		for cal in calibration:
			cal.resample(spectrum)

		# find best XCorr and amplitude adjustment
		best = None
		for cal in calibration:
			# best horizontal pixel shift
			shift = Xcorr( spectrum, cal, lag=lag)
			# build matrices with identical rows (len=npix-shift)
			if shift < 0:
				calmatrix = np.tile(cal.data[:shift], (trials, 1))
				objmatrix = np.tile(spectrum.data[-shift:], (trials, 1))
			elif shift > 0:
				calmatrix = np.tile(cal.data[shift:], (trials, 1))
				objmatrix = np.tile(spectrum.data[:-shift], (trials, 1))
			else:
				calmatrix = np.tile(cal.data, (trials, 1))
				objmatrix = np.tile(spectrum.data, (trials, 1))
			# amplitude matrix has identical columns
			size = np.shape(calmatrix)[1]
			ampmatrix = np.tile(amp, (size,1)).T

			# remove units for dimensionless operations
			calmatrix = calmatrix.value
			objmatrix = objmatrix.value

			# flip arrays for amplification
			diff = objmatrix - (1 - (1 - calmatrix) * ampmatrix)
			# compute the RMS for each trial
			rmsvector = np.sqrt(( diff**2 ).sum(axis=1) / size)
			if not best:
				# if first pass, assign to `best`
				best = ( rmsvector.min(), rmsvector.argmin(), shift, cal )
			elif rmsvector.min() < best[0]:
				# this fit is better, save it to `best`
				best = ( rmsvector.min(), rmsvector.argmin(), shift, cal )

		# results of calibration fitting
		index = best[1] # amplitude
		shift = best[2] # XCorr
		cal   = best[3] # which calibration spectrum

        # we can't update an attribute...
		update = spectrum.data.value
        
		# divide spectrum
		if shift < 0:
			update[-shift:] /= 1 - (1 - cal.data[:shift].value) * amp[index]
		elif shift > 0:
			update[:-shift] /= 1 - (1 - cal.data[shift:].value) * amp[index]
		else:
			update /= 1 - (1 - cal.data.value) * amp[index]

        # re-incorporate units
		spectrum.data = update * u.Unit(spectrum.yunits)

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise TelluricError('Inappropriate keyword arguments in '
		'Telluric.Correct().')
