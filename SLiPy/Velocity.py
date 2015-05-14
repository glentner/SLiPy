# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/SLiPy/Velocity.Py
"""
Radial velocity corrections for 1D spectra.
"""

from astropy.io import fits as pyfits
from astropy.constants import c

from .. import SlipyError
from ..astrolibpy.astrolib.helcorr import helcorr
from .Fits import Find, RFind
from .Observatory import Observatory
from .DataType import Spectrum
from ..Framework.Options import Options, OptionsError
from ..Framework.Display import Monitor, DisplayError

class VelocityError(SlipyError):
	"""
	Exception specific to the Velocity module
	"""
	pass

def HeaderInfo( fpath ):
	"""
	HeaderInfo( fpath ):

	Helper function of IrafInput.

	Return a formatted string containing the year, month, and day of
	observation from the FITS file with name 'fpath', as well as the
	univeral time of observation and the right ascension and declination
	of the target.
	"""
	try:
		with pyfits.open(fpath) as hdulist:
			ra       = ':'.join( hdulist[0].header['alpha'].split() )
			dec      = ':'.join( hdulist[0].header['delta'].split() )
			date_obs = hdulist[0].header['date-obs']
			date, UT = date_obs.split('T')
			year, month, day = [ x.lstrip('0') for x in date.split('-') ]

			info = '{:>5} {:>2} {:>2} {:>8} {:>11} {:>12}\n'.format(
				year, month, day, UT, ra, dec )

	except IOError as error:
		raise VelocityError('Failure in `{}` from HeaderInfo()'.format(fpath))

	return info

def IrafInput( *files, **kwargs):
	"""
	IrafInput( *args, **kwargs ):

	Build an input file for IRAF's rvcorrect task.

	`files` should be a list of FITS file names to build the output table for.
	The user can optionally specify a 'toplevel' directory to search
	(recursively!) under fitting the 'pattern' (default=*.fits). This results
	of this pattern search will be added to the list of file names in 'args'
	(if any given).

	kwargs = {
		'toplevel' : ''      , # search `toplevel` directory for files
		'pattern'  : '*.fits', # files under `toplevel` fitting `pattern`
		'recursive': False   , # search recusively under `toplevel`
		'outfile'  : ''        # write lines to file named `outfile`
	}
	"""
	try:
		# dictionary of options
		options = Options( kwargs,
			{
			'toplevel' : ''      , # search `toplevel` directory for files
			'pattern'  : '*.fits', # files under `toplevel` fitting `pattern`
			'recursive': False   , # search recursively under `toplevel`
			'outfile'  : ''        # write lines to file named `outfile`
			})

		# convert options types
		toplevel  = options('toplevel')
		pattern   = options('pattern')
		recursive = options('recursive')
		outfile   = options('outfile')

		files = list(files)

		if toplevel:
			# search for files matching `pattern`
			find   = RFind if recursive else Find
			files += find( toplevel, pattern )

		# get info from files
		info = [ HeaderInfo(fpath) for fpath in files ]

		if outfile:
			with open( outfile, 'w' ) as fp:
				fp.writelines( info )

		return info

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise VelocityError('Failed to construct table information in '
		'IrafInput().')

def HelioCorrect( obs, *spectra, **kwargs ):
	"""
	Perform heliocentric velocity corrects on `spectra` based on
	`obs`ervatory information (longitude, latitude, altitude) and the
	member attributes, ra (right ascension), dec (declination), and jd
	(julian date) from the `spectra`.
	"""
	try:

		# define function parameters
		options = Options( kwargs,
			{
				'verbose': False # display messages, progress
			})

		# assign parameters
		verbose = options('verbose')

		# check `obs` type
		if not issubclass( type(obs), Observatory):
			raise VelocityError('HelioCorrect() expects its first argument to '
			'be derived from the Observatory class.')
		elif ( not hasattr(obs, 'latitude') or not hasattr(obs,'longitude') or
			not hasattr(obs, 'altitude') ):
			raise VelocityError('HelioCorrect expects `obs`ervatory to have '
			'all three: latitude, longitude, and altitude attributes.')

		# check `spectra` arguments
		for a, spectrum in enumerate(spectra):
			if type(spectrum) is not Spectrum:
				raise VelocityError('HelioCorrect() expects all `spectrum` '
				'arguments to be of type Spectrum.')
			if not spectrum.ra or not spectrum.dec or not spectrum.jd:
				raise VelocityError('Spectrum {} lacks one or all of `ra`, '
				'`dec`, and `jd`; from HelioCorrect().'.format(a))

		if verbose:
				display = Monitor()
				print(' Running HelioCorrect on {} spectra ...'
				.format(len(spectra)))

		for a, spectrum in enumerate(spectra):
			# heliocentric velocity correction in km s^-1
			hcorr = helcorr(obs.longitude, obs.latitude, obs.altitude,
				spectrum.ra, spectrum.dec, spectrum.jd)[1]
			# apply correction to wave vector
			spectrum.wave -= spectrum.wave * 1000 * hcorr / c.value
			# show progress if desired
			if verbose: display.progress(a, len(spectra))

		# finalize progress bar (erase)
		if verbose: display.complete()

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise VelocityError('Failed to perform HelioCorrect() task.')

	except DisplayError as err:
		print(' --> DisplayError:', err)
		raise VelocityError('Exception from Display.Monitor in HelioCorrect().')


def BaryCorrect( obs, *spectra, **kwargs ):
	"""
	Perform barycentric velocity corrects on `spectra` based on
	`obs`ervatory information (longitude, latitude, altitude) and the
	member attributes, ra (right ascension), dec (declination), and jd
	(julian date) from the `spectra`.
	"""
	try:

		# define function parameters
		options = Options( kwargs,
			{
				'verbose': False # display messages, progress
			})

		# assign parameters
		verbose = options('verbose')

		# check `obs` type
		if not issubclass( type(obs), Observatory):
			raise VelocityError('HelioCorrect() expects its first argument to '
			'be derived from the Observatory class.')
		elif ( not hasattr(obs, 'latitude') or not hasattr(obs,'longitude') or
			not hasattr(obs, 'altitude') ):
			raise VelocityError('HelioCorrect expects `obs`ervatory to have '
			'all three: latitude, longitude, and altitude attributes.')

		# check `spectra` arguments
		for a, spectrum in enumerate(spectra):
			if type(spectrum) is not Spectrum:
				raise VelocityError('HelioCorrect() expects all `spectrum` '
				'arguments to be of type Spectrum.')
			if not spectrum.ra or not spectrum.dec or not spectrum.jd:
				raise VelocityError('Spectrum {} lacks one or all of `ra`, '
				'`dec`, and `jd`; from HelioCorrect().'.format(a))

		if verbose:
				display = Monitor()
				print(' Running HelioCorrect on {} spectra ...'
				.format(len(spectra)))

		for a, spectrum in enumerate(spectra):
			# heliocentric velocity correction in km s^-1
			hcorr = helcorr(obs.longitude, obs.latitude, obs.altitude,
				spectrum.ra, spectrum.dec, spectrum.jd)[0]
			# apply correction to wave vector
			spectrum.wave -= spectrum.wave * 1000 * hcorr / c.value
			# show progress if desired
			if verbose: display.progress(a, len(spectra))

		# finalize progress bar (erase)
		if verbose: display.complete()

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise VelocityError('Failed to perform HelioCorrect() task.')

	except DisplayError as err:
		print(' --> DisplayError:', err)
		raise VelocityError('Exception from Display.Monitor in HelioCorrect().')
