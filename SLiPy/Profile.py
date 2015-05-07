# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
"""
Profile fitting tasks for spectra.
"""
import numpy as np
from scipy.optimize import curve_fit

from .DataType import Spectrum, DataError
from ..Framework.Options import Options, OptionsError

class ProfileError(Exception):
	"""
	Exception specific to Profile module.
	"""
	pass

def igauss(x, *p):
	"""
	igauss(x, *p):

	Inverted Gaussian function -> y = 1 - A * exp( -(x-mu)**2 / (2*sigma**2) ).
	where `p` should be p = [A, mu, sigma].
	"""
	A, mu, sigma = p
	return 1 - np.exp( -(x - mu)**2 / (2 * sigma**2 ) )

def Fit(spectrum, xmin, xmax, **kwargs):
	"""
	Fit(spectrum, xmin, xmax, **kwargs):

	Given a `spectrum` of type Spectrum, and the `xmin`, `xmax` range specifying 
	the domain of the spectrum pertaining to the profile in question, fit a 
	Gaussian curve to the profile and return a result as type Spectrum.
	"""
	try:

		# function parameter defaults
		options = Options( kwargs,
			{
				'pad': 0 # expantion of domain for evaluation of the profile function.
			})

		# function parameter assignments 
		pad = options('pad')

		# check input arguments 
		if type(spectrum) is not Spectrum:
			raise ProfileError('Fit expects type Spectrum.')
		if xmin < spectrum.wave[0] or xmax > spectrum.wave[-1]:
			raise ProfileError('Out of domain for `spectrum`.')

		# make local copy of spectrum 
		local = spectrum.copy()

		# extract relavent data from spectrum 
		wave, data = local.wave, local.data 
		data = data[ wave[ wave < xmax ] > xmin ]
		wave = wave[ wave[ wave < xmax ] > xmin ]

		# initial guess of parameters 
		p0 = [ 1 - data.min(), wave[ data.argmin() ], data.std() ]

		# fit curve 
		coeff, var_matrix = curve_fit(igauss, wave, data, p0=p0)

		# create Spectrum of profile 
		x = local.wave 
		x = x[ x[ x < xmax + pad ] < xmin - pad ]
		profile_fit = Spectrum( igauss(x, *coeff) )
		profile_fit.wave = x

		return profile_fit

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise ProfileError('Failed keyword assignment in Fit().')

	except DataError as err:
		print(' --> DataError:', err)
		raise ProfileError('DataError in Fit().')

def Width(profile, **kwargs):
	"""
	"""
	pass
