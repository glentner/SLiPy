#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE file (GPLv2)

"""
usage: rvcorr.py @function [files] [**kwargs]

Radial velocity corrections for 1D spectra.
"""

import pyfits, os, sys, fnmatch

class InputError(Exception):
	"""
	Exception for bad arguments.
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

			info = '{:>5} {:>2} {:>2} {:>8} {:>11} {:>12}'.format( 
				year, month, day, UT, ra, dec )
		
	except IOError as error: 
		print( error )
		info = 'ERROR in {}'.format(fpath)
		print( info )
		pass
	
	return info

def IrafInput(*args, **kwargs):
	"""
	IrafInput( *args, **kwargs ):

	Build an input file for IRAF's rvcorrect task.

	'args' should be a list of FITS file names to build the output table for.
	The user can optionally specify a 'toplevel' directory to search 
	(recursively!) under fitting the 'pattern' (default=*.fits). This results
	of this pattern search will be added to the list of file names in 'args'
	(if any given). 

	if __name__ == '__main__', the function prints the contents of what would
	constitute 'outfile'; though the file is not actually created unless a 
	name is provided (e.g., outfile=rvcorr.dat). The printing can be suppressed
	via verbose=False.

	Options = { 
		'toplevel' : '', 'pattern' : '*.fits', 'outfile' : '', 'main' : False
	}
	"""
	# dictionary of options
	options = {
		'verbose'  : 'True',
		'toplevel' : '',
		'pattern'  : '*.fits',
		'outfile'  : '',
		'main'     : 'False'
	}

	for keyword in kwargs:
		# check for recognized options
		if keyword not in options:
			raise InputError(
				'{} is not a recognized keyword argument for IrafInput!'
				.format( keyword )
				)
		# assign to options
		options[keyword] = kwargs[keyword]

	# convert options types
	toplevel = options['toplevel']
	pattern  = options['pattern']
	outfile  = options['outfile']
	main     = options['main']
	if options['verbose'] == 'True':
		verbose = True
	elif options['verbose'] == 'False':
		verbose = False
	else: 
		raise InputError('verbose expects True or False')

	if toplevel:
		# match 'pattern' below 'toplevel'
		files = [
			os.path.join(toplevel, basename)
			for basename in fnmatch.filter( os.listdir(toplevel), pattern )
		]
	
	else: 
		files = [ ] # empty

	# all files to be searched
	files += list(args)
	
	# get info from files
	info = [ HeaderInfo(fpath) for fpath in files ]

	if outfile:
		with open( outfile, 'w' ) as fp:
			fp.writelines( info )

	if main and verbose:	
		for line in info:
			print(line)
	
	elif not main:
		return info

def Parse( *clargs ):
	"""
	Parse( *clargs ):

	Parse clargs = sys.argv to seperate args and kwargs.
	Returns args, kwargs.
	"""
	# args should not have an assignment
	args = [ x for x in clargs if '=' not in x ]

	# build dictionary of remaining kwargs
	kwargs = { 
		key:value for key, value in [ 
			arg.split('=') for arg in set(clargs) - set(args) 
		] 
	}
	
	return args, kwargs


if __name__ == '__main__':

	if len( sys.argv ) < 2:
		# print usage
		print( __doc__ )
		sys.exit(0)

	# executable functions from command line
	executable = {
		'IrafInput' : IrafInput.__doc__ 
	}

	# retrieve name of executable
	functioncall = sys.argv[1].split('@')
	if len(functioncall) != 2 or functioncall[0] != '':
		print(
			'\n<{}> was not understood as a valid function call. '
			' Execute {} with out arguments for usage details.\n'
			.format( sys.argv[1], os.path.basename(sys.argv[0]) )
		)
		sys.exit(1)
	
	else: function = functioncall[1]
	
	# check for valid function call
	if function not in executable:
		print(
			'\n<{}> was not understood as a valid function call. '
			'Execute {} with out arguments for usage details.\n'
			.format( function, os.path.basename(sys.argv[0]) )
		)
		sys.exit(1)

	if len(sys.argv) < 3 or sys.argv[2] == 'help':
		# show function usage
		print( executable[function] )
		sys.exit(0)

	try:
		# select out arguments and keywords
		args, kwargs = Parse( *sys.argv[2:] )
	
	except ValueError as error:
		print( error )
		print( '--> InputError: syntax error in your key word assignments!' )
		sys.exit(1)

	# choose function
	if function == 'IrafInput':
		EXE = IrafInput

	# Run script
	EXE( *args, **kwargs )
