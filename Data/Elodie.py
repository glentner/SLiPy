# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Data/Elodie.py
"""
Methods for data retrieval from the Elodie Archive.
"""

import os, numpy as np
from urllib.request import urlopen
from ..Framework.Options import Options, OptionsError 

class ElodieError(Exception):
	"""
	Exception specific to Elodie module.
	"""
	pass

class Archive:
	"""
	Import and parse ascii catalog of Elodie archive files. The complete
	archive is stored in the member `data`. It's organized in a dictionary
	by unique target names. The reduced archive of contains identified `HD`, 
	`BD`, `GC`, and `GJ` objects, choosing the file pertaining to the spectra
	with the highest signal-to-noise ratio available.
	"""
	def __init__(self, **kwargs):
		try:

			# default input file 
			default_infile = os.path.join( os.path.dirname(__file__), 
				'archives/elodie.csv')

			# parameter defaults
			options = Options( kwargs, 
				{
					'infile' : default_infile       , # path to input file 
					'catalog': ['HD','BD','GC','GJ']  # catalogues to keep
				})

			# parameter assignments 
			infile  = options('infile')
			catalog = options('catalog')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise ElodieError('Failed keyword assignment from Archive.__init__')

		# import data from archive
		with open(infile, 'r') as archive:
			data = [ line.split(',') for line in archive.readlines() ]

		# strip elements of white space 
		data = [ [x.strip() for x in line] for line in data ]

		# build empty dictionary of unique names with empty lists
		targets = { name:[] for name in set([ line[0] for line in data ]) }

		# compile list of spectra files organized by identifier w/ S/N
		for line in data:
			targets[line[0]].append( line[2:] )
		
		# reject files from spectra not in catalogues
		targets = { k:v for k, v in targets.items() if k[:2] in catalog }

		files = {}
		for target, options in targets.items():
			# choose best S/N
			idx = np.array( np.array(options)[:,-1], dtype=np.int_ ).argmax()
			files[target] = options[idx][0]
			
		# members 
		self.data  = data 
		self.files = files 

def Script(filename, pipeline=''):
	"""
	Construct url script for Elodie archive given `filename` and optionally
	`pipeline` instructions (e.g., `&z=wrs|fca[1,nor]`).
	"""
	return ''.join(['http://atlas.obs-hp.fr/elodie/E.cgi?&c=i&o=',
		filename, pipeline, '&a=mime:application/x-fits'])

def Download( *files, **kwargs ):
	"""
	Download `files` from Elodie archive via url scripts. The spectra can be 
	further reduced via Elodie`s pipeline with the following options.

	kwargs = {
			'verbose'  : True           , # display messages, progress
			'resample' : (min, max, res), # resample spectra (no default)
			'normalize': True           , # continuum normalization
			'outpath'  : './'           , # directory for downloaded files
			'names'    : []               # alternative output names for `files`
		}
	"""
	try:

		# function parameter defaults 
		options = Options( kwargs,
			{
				'verbose'  : True   , # display messages, progress
				'resample' : (0,0,0), # no default 
				'normalize': True   , # continuum normalization 
				'outpath'  : './'   , # directory for downloaded files
				'names'    : []       # alternative output names for `files`
			})

		# function parameter assignments 
		verbose   = options('verbose')
		resample  = options('resample')
		normalize = options('normalize')
		outpath   = options('outpath')

		# check for `resampled` assignment 
		if 'resample' not in kwargs:
			resample = None

		# set default names 
		if not names:
			names = '-'.join(fname.split(':'))

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise ElodieError('Failed keyword assignment from Download().')

	if not resample and not normalize:
		pipeline = ''
	else:
		pipeline = '&z=wrs'

	if normalize:
		pipeline += '|fca[1,nor]'

	if resample:
		pipeline += '|wrs[1' + ','.join(resample) + ']'

	if verbose:
		display = Monitor()
		nfiles  = len(files)
		print('\n Downloading {} files from Elodie ...'.format(nfiles))

	for a, spectra in enumerate(files):

		# download file
		with urlopen( Script(spectra) ) as response, open(
			os.path.join(outpath, names[a]), 'wb') as outfile:
			    shutil.copyfileobj(response, outfile)

		# show progress 
		if verbose: display.progress(a, nfiles)

	if verbose:
		display.complete()

