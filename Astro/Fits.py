#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/Astro/Fits.py 
"""
Fits:

Module for importing data and header information from FITS image files.
"""
import os, sys, pyfits, fnmatch
from Python import BaseError
from Python.General.Interface import Parse, CommandError
from Python.General.Options import Options, OptionsError
from Python.General.Display import Display, DisplayError
from Python.Astro.DataType import Spectrum, DataError
from Python.Astro.Simbad import Position, Distance, \
		Sptype, IDList, SimbadError

class FitsError(BaseError):
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
		raise FitsError('`{}` does not name a directory.'
				.format(toplevel))

	return [
			os.path.join(toplevel, filename)
			for filename in fnmatch.filter(os.listdir(toplevel), pattern)
		]

def RFind(toplevel, pattern):
	"""
	RFind(toplevel, pattern):

	Recursively search for paths below `toplevel` fitting `pattern`.
	"""
	if not os.path.isdir(toplevel):
		raise FitsError('`{}` does not name a directory.'
				.format(toplevel))

	return [
			os.path.join(dirpath, filename)
			for dirpath, dirnames, filenames in os.walk(toplevel)
			for filename in fnmatch.filter(filenames, pattern)
		]

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
			print(' Importing data from {} Fits files ...'
					.format(nfiles) )
			for a, filename in enumerate(files):
				display.progress(a, nfiles)
				data.append(
						Spectrum(filename, wavecal=wavecal, crpix1=crpix1,
							crval1=crval1, cdelt1=cdelt1)	
					)
			display.complete()
			return data

		# import spectra 'silently'
		return [
				Spectrum(filename, wavecal=wavecal, crpix1=crpix1,
					crval1=crval1, cdelt1=cdelt1) for filename in files
			]
	
	except OptionsError as err:
		print(' --> OptionsError:', err.msg)
		raise FitsError('Data retrieval failure.')

	except DataError as err:
		print(' --> DataError:', err.msg)
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
		print(' --> OptionsError:', err.msg)
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
			print(' Reading object names for {} Fits files ...'
					.format(nfiles))
			obj_ids = []
			for a, name in enumerate(files):
				display.progress(a, nfiles)
				obj_ids.append( Header(name, 'object') )
		
			display.complete()

		else: obj_ids = [
					Header(name, 'object') for name in files
				]

		if verbose:
			# query for `attribute` iteratively
			print(' Searching for `{}`s with SIMBAD ...'
					.format(attribute))
			results = []
			for a, obj in enumerate(obj_ids):
				display.progress(a, nfiles)
				results.append(
						SimbadSearch[attribute](obj)	
					)
			display.complete()

		else: results = [
					SimbadSearch[attribute](obj)
					for obj in obj_ids
				]

		if is_main:
			formatted = {
					'Position':'{0:.2f} {1:.2f}',
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
		print(' --> OptionsError:', err.msg)
		raise FitsError('Failed assignment for Search().')

	except SimbadError as err:
		print(' --> SimbadError:', err.msg)
		raise FitsError('Simbad failed.')

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
		print(' --> CommandError:', err.msg)
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
