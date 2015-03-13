#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/Fits.py 
"""
Fits - FITS file handling module.
"""
import os, sys, fnmatch
from astropy.io import fits as pyfits 
from numbers import Number
from ..Framework.Command import Parse, CommandError
from ..Framework.Options import Options, OptionsError
from ..Framework.Display import Monitor, DisplayError
from .DataType import Spectrum, DataError
from .Simbad import Position, Distance, Sptype, IDList, SimbadError

class FitsError(Exception):
	"""
	Exception for Fits module.
	"""
	pass

def Find(toplevel, pattern):
	"""
	Find(toplevel, pattern):

	Search for file paths below `toplevel` fitting `pattern`.
	"""
	if not os.path.isdir(toplevel):
		raise FitsError('`{}` does not name a directory.'.format(toplevel))

	return [ os.path.join(toplevel, filename)
		for filename in fnmatch.filter(os.listdir(toplevel), pattern) ]

def RFind(toplevel, pattern):
	"""
	RFind(toplevel, pattern):

	Recursively search for paths below `toplevel` fitting `pattern`.
	"""
	if not os.path.isdir(toplevel):
		raise FitsError('`{}` does not name a directory.'.format(toplevel))

	return [ os.path.join(dirpath, filename)
		for dirpath, dirnames, filenames in os.walk(toplevel)
		for filename in fnmatch.filter(filenames, pattern) ]

def GetData( *files, **kwargs ):
	"""
	GetData( *files, **kwargs ):

	Import data from FITS `files`.

	kwargs = {
			verbose   : True    , # display messages, progress
			toplevel  : ''      , # request import from directory `toplevel`
			pattern   : '*.fits', # pattern matching with `toplevel`
			recursive : False   , # search recursively below `toplevel`
			wavecal   : True    , # fit wavelength vector to data 
			crpix1    : 'crpix1', # reference pixel header keyword 
			crval1    : 'crval1', # value at reference pixel 
			cdelt1    : 'cdelt1', # resolution (delta lambda)
		}
	"""
	try:
		# convert `files` to list 
		files = list(files)

		# available key word arguments 
		options = Options( kwargs,
			{
				'verbose'  : True    , # display messages, progress
				'toplevel' : ''      , # request import from `toplevel` dir
				'pattern'  : '*.fits', # pattern matching with `toplevel`
				'recursive': False   , # search recursively below `toplevel`
				'wavecal'  : True    , # fit wavelength vector to data 
				'crpix1'   : 'crpix1', # reference pixel header keyword 
				'crval1'   : 'crval1', # value at reference pixel 
				'cdelt1'   : 'cdelt1', # resolution (delta lambda)
			})
		
		# assignment options 
		verbose   = options('verbose')
		toplevel  = options('toplevel')
		pattern   = options('pattern')
		recursive = options('recursive')
		wavecal   = options('wavecal')
		crpix1    = options('crpix1')
		crval1    = options('crval1')
		cdelt1    = options('cdelt1')
		
		if toplevel:
			# search for files matching `pattern`
			find   = RFind if recursive else Find 
			files += find( toplevel, pattern )
		
		if verbose:
			# import iteratively, displaying progress
			display = Monitor()
			nfiles  = len(files)
			data    = []
			print(' Importing data from {} Fits files ...'.format(nfiles) )
			for a, filename in enumerate(files):
				display.progress(a, nfiles)
				data.append( Spectrum(filename, wavecal=wavecal,
					crpix1=crpix1, crval1=crval1, cdelt1=cdelt1) )

			display.complete()
			return data

		# import spectra 'silently'
		return [ Spectrum(filename, wavecal=wavecal, crpix1=crpix1,
			crval1=crval1, cdelt1=cdelt1) for filename in files ]
	
	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise FitsError('Data retrieval failure.')

	except DataError as err:
		print(' --> DataError:', err)
		raise FitsError('Failed to construct spectrum.')

def Header( filename, keyword, **kwargs ):
	"""
	Header( filename, keyword, **kwargs ):

	Retrieve `keyword` from Fits `filename`.
	"""
	try:

		options = Options( kwargs,
			{
				'is_main' : False # reserved for calls from Main()
			})

		is_main = options('is_main')
	
		with pyfits.open(filename) as hdulist:
			element = hdulist[0].header[keyword]

		if is_main:
			print( element )
			return

		else: return element

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise FitsError('Failed keyword assignment in Header().')

	except KeyError as key:
		raise FitsError('Header element `{}` was not accessible '
			'from `{}`'.format(keyword, filename))
	
def Search( *files, **kwargs ):
	"""
	Search( *files, **kwargs ):

	Exctract object names from Fits `files` and use Simbad.py 
	to resolve the `attribute` (a required keyword argument)
	from the SIMBAD astronomical database.

	kwargs = {
			verbose   : True    , # display messages, progress 
			toplevel  : ''      , # search under `toplevel` directory 
			pattern   : '*.fits', # for files under `toplevel`
			recursive : False   , # search recusively under `toplevel` 
			attribute : ''      , # attribute to search for (no default)
		}
	"""
	try:
		
		# convert `files` to list 
		files = list(files)

		# available keyword arguments
		options = Options( kwargs,
			{
				'verbose'   : True    , # display messages, progress 
				'toplevel'  : ''      , # search under `toplevel` directory 
				'pattern'   : '*.fits', # for files under `toplevel`
				'recursive' : False   , # search recusively under `toplevel` 
				'attribute' : ''      , # attribute to search for (no default)
				'is_main'   : False     # reserved for calls from Main()
			})

		# assign parameters 
		verbose   = options('verbose')
		toplevel  = options('toplevel')
		pattern   = options('pattern')
		recursive = options('recursive')
		attribute = options('attribute')
		is_main   = options('is_main')

		# Available search functions from Simbad.py
		SimbadSearch = {
				'Position': Position, # ra, dec	(degrees)
				'Distance': Distance, # in parsecs
				'Sptype'  : Sptype  , # spectral type
				'IDList'  : IDList    # alternate IDs
			}

		if not attribute:
			raise FitsError('An `attribute` must be specified for Search().')
		
		if attribute not in SimbadSearch:
			raise FitsError('`{}` is not an available search criteria.'
					.format(attribute))

		if toplevel:
			# search for files in `toplevel` directory
			find   = RFind if recursive else Find
			files += find(toplevel, pattern)

		nfiles = len(files)
		display = Monitor()
		
		if verbose:
			# read object names iteratively
			print(' Reading object names for {} Fits files ...'.format(nfiles))
			obj_ids = []
			for a, name in enumerate(files):
				display.progress(a, nfiles)
				obj_ids.append( Header(name, 'object') )
		
			display.complete()

		else: obj_ids = [ Header(name, 'object') for name in files ]

		if verbose:
			# query for `attribute` iteratively
			print(' Searching for `{}`s with SIMBAD ...'.format(attribute))
			results = []
			for a, obj in enumerate(obj_ids):
				display.progress(a, nfiles)
				results.append( SimbadSearch[attribute](obj) )
			
			display.complete()

		else: results = [ SimbadSearch[attribute](obj) for obj in obj_ids ]

		if is_main:
			formatted = {
					'Position':'{1:.2f} {1:.2f}',
					'Distance':'{0:.2f}',
					'Sptype'  : '{}'
				}
			for item in results:
				if type(item) is list:
					print( formatted[attribute].format(*item))
				else:
					print( formatted[attribute].format(item))

		else: return results

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise FitsError('Failed assignment for Search().')

	except SimbadError as err:
		print(' --> SimbadError:', err)
		raise FitsError('Simbad failed.')

def PositionSort( center, radius, *files, **kwargs ):
	"""
	PositionSort( *files, **kwargs ):

	Return a list of files from `files` that lie in a `radius` (in degrees)
	from `center`, based on the `ra` (right ascension) and `dec` (declination).
	
	kwargs = {
			'ra'       : 'pos1'  , # header element for right ascension 
			'dec'      : 'pos2'  , # header element for declination
			'obj'      : 'object', # header element for object id
			'raconvert': True    , # convert decimal hours to decimal degrees
			'verbose'  : True    , # display messages, progress
			'toplevel' : ''      , # `toplevel` directory to look for files in 
			'recursive': False   , # search `recursive`ly below `toplevel`
			'pattern'  : '*.fits', # glob `pattern` for file search 
			'useSimbad': False     # use Simbad instead of header elements
		}
	"""
	try:
		# function parameter defaults
		options = Options( kwargs,
			{
				'ra'       : 'pos1'  , # header element for right ascension 
				'dec'      : 'pos2'  , # header element for declination
				'obj'      : 'object', # header element for object id
				'raconvert': True    , # convert decimal hours to decimal degrees
				'verbose'  : True    , # display messages, progress
				'toplevel' : ''      , # `toplevel` directory for file search
				'recursive': False   , # search `recursive`ly below `toplevel`
				'pattern'  : '*.fits', # glob `pattern` for file search 
				'useSimbad': False     # use Simbad instead of header elements
			})

		# function parameter assignments
		ra        = options('ra')
		dec       = options('dec')
		obj       = options('obj')
		raconvert = options('raconvert')
		verbose   = options('verbose')
		toplevel  = options('toplevel')
		recursive = options('recursive')
		pattern   = options('pattern')
		useSimbad = options('useSimbad')

		# check arguments 
		if not hasattr( center, '__iter__'):
			raise FitsError('PositionSort() expects `center` argument to be '
			'`iterable` and have two elements.')
		if len(center) != 2:
			raise FitsError('PositionSort() expects `center` argument to have '
			'exactly two elements.')
		if not isinstance( radius, Number ):
			raise FitsError('PositionSort() expects `radius` argument to '
			'be a `Number`.')
		for a, f in enumerate(files):
			if not isinstance(f, str):
				raise FitsError('PositionSort() expects `str` like arguments '
				'for all `files` (from argument {})'.format(a))

		# convert `files` to list type 
		files = list(files)

		# look under `toplevel` if requested
		if toplevel:
			find = RFind if recursive else Find
			files += find(toplevel, pattern)

		if verbose:
			# create display object 
			display = Monitor()
			nfiles  = len(files)

		# initialize blank lists
		pos1, pos2, = [], []

		if not useSimbad:	
			if verbose: print(' Retrieving {} positions from files ... '
				.format(nfiles))
			# check file headers for requested information
			for a, fitsfile in enumerate(files):
				pos1.append( Header(fitsfile, ra)  )
				pos2.append( Header(fitsfile, dec) )
				if verbose: display.progress( a, nfiles )

		else:
			# use the Simbad module to search for positions 
			if verbose: print(' Retrieving {} positions from SIMBAD ... '
				.format(nfiles))
			for a, fitsfile in enumerate(files):
				pos = Position( Header(fitsfile, obj) )
				pos1.append( pos[0] )
				pos2.append( pos[1] )
				if verbose: display.progress(a, nfiles)

		# erase progress bar
		if verbose: 
			display.complete()
			print(' Compiling list of files ... ')

		# keep files for targets within range
		keepers = [ f for p1, p2, f in zip(pos1, pos2, files) 
			if abs(p1 - center[0]) < radius and abs(p2 - center[1]) < radius ]

		# account for p1 ~ 0 && center ~ 360 like comparisons
		keepers += [ f for p1, p2, f in zip(pos1, pos2, files)
			if abs(p1 + 360 - center[0]) < radius and 
			abs(p2 - center[1]) < radius ]

		# account for p1 ~ 360  && center ~ 0 like comparisons
		keepers += [ f for p1, p2, f in zip(pos1, pos2, files)
			if abs(p1 - center[0] - 360) < radius and 
			abs(p2 - center[1]) < radius ]

		if verbose:
			print('\033[1A\r Compiling list of files ... done')

		# exclude any potential double countings
		return list( set(keepers) )


	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise FitsError('Failed keyword assignment in PositionSort().')

def Main( clargs ):
	"""
	Main function. See __doc__ for details.
	"""
	
	if len(clargs) < 2:
		# show usage 
		print( __doc__ )
		return 0
	
	# Available functions for execution
	executable = {
			'Header' : Header, # Header function
			'Search' : Search  # Search function
		}

	try:

		# Parse command line arguments 
		function, args, kwargs = Parse( clargs[1:] )
		if not args and not kwargs:
			# show function usage
			print( executable[function].__doc__ )
			return 0

		# run execution 
		executable[function]( *args, is_main=True, **kwargs )
		return 0

	except CommandError as err:
		print(' --> CommandError:', err)
		return 1

	except KeyError as key:
		print(' --> {} was not a recognized function.'.format(key))
		return 1

	except FitsError as err:
		# don't let uncaught self exception pass if from main.
		print(' --> FitsError:', err.msg)
		return 1

	except Exception as err:
		print(' --> Unrecognized error from Fits module.')
		print(' --> Exception: `{}`'.format(err))
		return 1


if __name__ == '__main__':
	# call Main function, exit 0 or 1
	sys.exit( Main(sys.argv) )
