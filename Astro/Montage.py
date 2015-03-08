# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/Montage.py 
"""
This module makes use of the `Montage` mosaic tools from caltech, see:
http://montage.ipac.caltech.edu/

The user should have Montage`s executables available on their path.
"""

import os, shutil as sh, numpy as np 
from subprocess import check_output as call, CalledProcessError
from sys import stdout
from numbers import Number
from ..Framework.Options import Options, OptionsError
from ..Framework.Display import Monitor, DisplayError


class MontageError(Exception):
	"""
	Exception specific to the Montage module
	"""
	pass

def SolveGrid( sides, grid ):
	"""
	SolveGrid( sides, grid ):

	Helper function for the Field and SubField classes. Both `sides` and `grid`
	need to be array-like and of length two. `sides` is the side length of the 
	field in decimal degrees in right ascension and declination respectively.
	`grid` specifies the subdivision along these axis (e.g., (2,2) says 2x2).
	
	The user should mindful of their choices. If the side lengths cannot be
	subdivided into well-behaved (rational) segments, higher decimal places
	will be lossed in the SubField.ArchiveList() task resulting in small 
	gaps in the mosaic.
	"""
	# check arguments 
	if not hasattr(sides, '__iter__') or not hasattr(grid, '__iter__'):
		raise MontageError('Grid() expects both arguments to be array-like.')
	if len(sides) != 2 or len(grid) != 2:
		raise MontageError('Grid() expects both arguments to have length two.')

	# grid `site` centers in the horizontal axis
	ra_left_site    = -sides[0] / 2 + 0.5 * sides[0] / grid[0]
	ra_right_site   =  sides[0] / 2 - 0.5 * sides[0] / grid[0]
	ra_site_centers = np.linspace( ra_left_site, ra_right_site, grid[0] )

	# grid `site` centers in the vertical axis 
	dec_bottom_site  = -sides[1] / 2 + 0.5 * sides[1] / grid[1]
	dec_top_site     =  sides[1] / 2 - 0.5 * sides[1] / grid[1]
	dec_site_centers = np.linspace( dec_bottom_site, dec_top_site, grid[1] )

	return ra_site_centers, dec_site_centers

def Mosaic(resolution, *folders, **kwargs):
	"""
	Mosaic(resolution, *folders, **kwargs):

	Conduct standard build procedures for all `folders`. `resolution` is the 
	number of pixels per degree for the output image. Note: `folders` should 
	be absolute paths.

	kwargs = {
			verbose : True, # display messages, progress
			bkmodel : True  # model and correct for background effects
		}
	"""
	try:
		# function parameter defaults
		options = Options( kwargs,
			{
				'verbose': True, # display messages, progress
				'bkmodel': True  # model andd correct for background effects 
			})

		# function parameter assignments 
		verbose = options('verbose')
		bkmodel = options('bkmodel')

		# build mosaics at all sites 
		for a, folder in enumerate(folders):
			
			# change directories 
			os.chdir(folder)

			# display message
			if verbose: 
				stdout.write('\n Building mosaic for {}: {} / {} ... \n'
					.format(os.path.basename(folder), a + 1, len(folders)))
				stdout.write( ' ' + '-' * 70 + '\n' )
				stdout.write(' Generating image meta-data table ... ')
				stdout.flush()
			
			# generate image meta-data table
			output = call(['mImgtbl','images','images.tbl']).decode('utf-8')
			if 'ERROR' in output: raise MontageError('Failed `mImgtbl` from `{}`'
				.format(folder))

			if verbose: 
				stdout.write('done.\n Generating FITS header template ... ')
				stdout.flush()

			# create mosaic FITS header template
			output = call(['mMakeHdr','-p','{}'.format(1 / resolution),
				'-n','images.tbl','template.hdr']).decode('utf-8')
			if 'ERROR' in output: raise MontageError('Failed `mMakeHdr` from '
				'`{}`'.format(folder))
			
			if verbose: 
				stdout.write('done\n Reprojecting images ... ')
				stdout.flush()

			# reproject images 
			output = call(['mProjExec','-p','images','images.tbl',
				'template.hdr','projected','stats.tbl']).decode('utf-8')
			if 'ERROR' in output: raise MontageError('Failed `mProjExec` in '
				'`{}`'.format(folder))

			if verbose: 
				stdout.write('done\n Generating new image meta-data table '
					'for projected images ... ')
				stdout.flush()

			# create new meta-data table for reprojected images 
			output = call(['mImgtbl','projected','proj-images.tbl'
				]).decode('utf-8')
			if 'ERROR' in output: raise MontageError('Failed `mImgtbl` in '
				'`{}`'.format(folder))
			
			if not bkmodel:
				# simply co-add images 
				if verbose: 
					stdout.write('done\n Co-adding images ... ')
					stdout.flush()
				output = call(['mAdd','-p','projected','proj-images.tbl',
					'template.hdr','final/mosaic.fits']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mAdd` in '
					'`{}`'.format(folder))

			else:
				# Fit overlaps for background corrections
				if verbose: 
					stdout.write('done\n Fitting overlaps for background '
						'corrections ... ')
					stdout.flush()
				output = call(['mOverlaps','proj-images.tbl','diffs.tbl'
					]).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mOverlaps` in '
					'`{}`'.format(folder))

				# perform background subtractions on overlaps
				if verbose: 
					stdout.write('done\n Performing background subtractions '
						'on overlaps ... ')
					stdout.flush()
				output = call(['mDiffExec','-p','projected','diffs.tbl',
					'template.hdr','differences']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mDiffExec` in '
					'`{}`'.format(folder))
	
				# computing plane-fitting coefficients
				if verbose: 
					stdout.write('done\n Computing plane-fitting coefficients ... ')
					stdout.flush()
				output = call(['mFitExec','diffs.tbl','fits.tbl',
					'differences']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mFitExec` in '
					'`{}`'.format(folder))

				# create table of background corrections
				if verbose: 
					stdout.write('done\n Creating table of background '
						'corrections ... ')
					stdout.flush()
				output = call(['mBgModel','proj-images.tbl','fits.tbl',
					'corrections.tbl']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mBgModel` in '
					'`{}`'.format(folder))

				# apply background matching to reprojected images 
				if verbose: 
					stdout.write('done\n Applying background matching to '
						'reprojected images ... ')
					stdout.flush()
				output = call(['mBgExec','-p','projected','proj-images.tbl',
					'corrections.tbl','corrected']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mBgExec` in '
					'`{}`'.format(folder))

				# co-add images for final mosaic 
				if verbose: 
					stdout.write('done\n Co-adding corrected images ... ')
					stdout.flush()
				output = call(['mAdd','-p','corrected','proj-images.tbl',
					'template.hdr','final/mosaic.fits']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mAdd` in '
					'`{}`'.format(folder))

			# finished mosaic
			if verbose:
				stdout.write('done\n')
				stdout.flush()

	except CalledProcessError as err:
		print(' --> CalledProcessError:', err)
		raise MontageError('Failed process from Mosaic().')

class SubField:
	"""
	SubField( center, sides, grid, **kwargs ):
	"""
	def __init__(self, center, sides, grid, **kwargs ):
		"""
		Create `site` grid for SubField.
		"""
		try:
			# function parameter defaults
			options = Options( kwargs,
				{
					'survey': 'DSS'  , # DSS, SDSS, 2MASS
					'band'  : 'DSS2B', # filter for `survey`, see `bands` dict
				})

			# function parameter assigments
			survey = options('survey').upper()
			band   = options('band').upper()

			# available bands for each survey 
			bands = {
					# Two-Micron All-Sky Survey
					'2MASS': ['J', 'H', 'K'],
					
					# Sloan Digital Sky Survey
					'SDSS': ['U','G','R','I','Z'],

					# STScI Digitized Sky Survey 
					'DSS': ['DSS1B','DSS1R','DSS2B','DSS2R','DSS2IR','Quick-V']
				}

			# check for appropriate survey, band 
			if survey not in bands:
				raise MontageError('`{}` was not a recognized survey for '
				'the Montage Python module.'.format(survey))
			if band not in bands[survey]:
				raise MontageError('`{}` was not a recognized filter band '
				'for the `{}` survey.'.format(band, survey))

			# check arguments 
			if ( not hasattr(center, '__iter__') or 
				not hasattr(sides, '__iter__') or not hasattr(grid, '__iter__') ):
				raise MontageError('SubField() expects array-like arguments for '
				'`center`, `sides`, and `grid` arguments.')
			if len(center) != 2 or len(sides) != 2 or len(grid) != 2:
				raise MontageError('SubField() expects `center`, `sides` and '
				'`grid` arguments to have length two.')

			# SolveGrid()
			ra_site_centers, dec_site_centers = SolveGrid(sides, grid)
			# relative to SubField `center`:
			ra_site_centers  += center[0]
			dec_site_centers += center[1] 

			# record number of `site`s along each axis
			self.num_ra_sites  = grid[0]
			self.num_dec_sites = grid[1]

			# build arguments for subprocess call
			self.archive_command_list = [
					['mArchiveList', survey, band, '{:.2f} {:.2f}'.format(ra_site, 
						dec_site), str(sides[0]/grid[0]), str(sides[1]/grid[1]), 
						'remote.tbl']
					for ra_site in ra_site_centers for dec_site in dec_site_centers
				]
		
			# make current directory `home`
			self.location = os.path.abspath('.')

			# new tree structure
			self.folders = [os.path.join(self.location,'Site-{:03d}'.format(a+1)) 
				for a in range(len(self.archive_command_list)) ]
		
			# initialize folder structure with `images` directory
			for folder in self.folders:
				abspath = os.path.join(folder,'images')
				if not os.path.exists(abspath):
					os.makedirs(abspath)

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in SubField().')

	def ArchiveList(self, **kwargs):
		"""
		Run the `mArchiveList` command on the `site` grid.
		"""
		try:
			# function parameter defaults
			options = Options( kwargs,
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments
			verbose = options('verbose')

			for a, command in enumerate(self.archive_command_list):
				
				# navigate to `images` directory 
				os.chdir( os.path.join(self.folders[a],'images') )
				
				if verbose:
					stdout.write('\n Running `mArchiveList` on {} ... '.format(
						os.path.basename(self.folders[a])))
					stdout.flush()
				
				# run `mArchiveList`
				output = call(command).decode('utf-8')
				if 'ERROR' in output or 'count="0"' in output:
					raise MontageError('Failed archive list from archive() '
					'(command: {}), (output: {}).'.format(command, output))

				if verbose: stdout.write('done')

			if verbose:
				stdout.write('\n')
				stdout.flush()

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigned from SubField.exec().')

		except CalledProcessError as err:
			print(' --> CalledProcessError:', err)
			raise MontageError('`mArchiveList` returned exit status 1.')

	def ArchiveExec(self, **kwargs):
		"""
		Run `mArchiveExec` on each `site` in the SubField.
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs,
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments
			verbose = options('verbose')

			for a, folder in enumerate(self.folders):
				
				# navigate to site folder 
				os.chdir( os.path.join(folder,'images') )
				
				if verbose:
					stdout.write('\n Running `mArchiveExec` on {} ... '.format(
						os.path.basename(folder)))
					stdout.flush()

				# run `mArchiveExec`
				output = call(['mArchiveExec','remote.tbl']).decode('utf-8')
				if 'ERROR' in output: raise MontageError('Failed `mArchiveExec` '
					'in folder `{}`. Output: {}'.format(folder, output))
				
				if verbose: stdout.write('done')

			if verbose: 
				stdout.write('\n')
				stdout.flush()

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigment in ArchiveExec().')

		except CalledProcessError as err:
			print(' --> CalledProcessError:', err)
			raise MontageError('`mArchiveExec` returned exit status 1.')

	def Build(self, resolution, **kwargs):
		"""
		Run the build process for the `sites` in this SubField. See the 
		Montage.Mosaic() function documentation.
		"""
		try:
			# function parameter options 
			options = Options( kwargs,
				{
					'verbose':True, # display message, progress 
					'bkmodel':True  # run background modelling procedure
				})
			
			# function parameter assignments 
			verbose = options('verbose')
			bkmodel = options('bkmodel')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in SubField.Build().')
	
		if verbose: 
			stdout.write(' Setting up folder structure ... ')
			stdout.flush()
		
		# setup folder structure 
		for subdir in ['corrected','projected','differences','final']:
			for folder in self.folders:
				abspath = os.path.join(folder,subdir)
				if not os.path.exists(abspath):
					os.makedirs(abspath)
		
		if verbose: 
			stdout.write('done\n')
			stdout.flush()

		# run Mosaic() on all `site`s
		Mosaic( resolution, *self.folders, verbose=verbose, bkmodel=bkmodel )

	def Collect(self, **kwargs):
		"""
		Collect( **kwargs ):

		Collect the mosaics from all `site` locations into a master `images`
		folder. The only `kwarg` is `verbose` (default: True).
		"""
		try:
			
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True # display messages, progress 
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in SubField.Collect().')

		# change to SubField `location` directory 
		os.chdir(self.location)

		# check that we have a finished mosaic at each `site`
		for folder in self.folders:
			if not os.path.exists( os.path.join(folder, 'final/mosaic.fits') ):
				raise MontageError('No `mosaic.fits` file in `{}`.'.format(folder))

		# ensure path to master/images/
		master_images = os.path.join(self.location, 'master/images')
		if not os.path.exists(master_images):
			os.makedirs(master_images)

		for a, folder in enumerate(self.folders):

			if verbose:
				stdout.write('\n Copying image from {} ... '.format(
					os.path.basename(folder)))
				stdout.flush()

			# Copy over image 
			sh.copy( os.path.join(folder,'final/mosaic.fits'), os.path.join(
				master_images, 'mosaic-{}.fits'.format(a + 1)) )

			if verbose: stdout.write('done')

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def Merge(self, resolution, **kwargs):
		"""
		Merge(resolution, **kwargs ):

		Merge all `site` mosaics into a single master SubField mosaic. The only 
		keyword options are `verbose` (default: True) and `bkmodel` (default:
		True). See Montage.Mosaic().
		"""
		try:
			
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True, # display messages, progress 
					'bkmodel': True  # model and correct for background effects
				})

			# function parameter assignments 
			verbose = options('verbose') 
			bkmodel = options('bkmodel')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigment in SubField.Merge().')

		# check for master directory
		master_dir = os.path.join(self.location,'master')
		if not os.path.exists(master_dir):
			raise MontageError('No `master` directory detected for SubField '
				'`{}`.'.format(self.location))
		
		# change directories 
		os.chdir(master_dir)

		# create folder structure 
		for subdir in ['corrected','projected','differences','final']:
			path = os.path.join(master_dir,subdir)
			if not os.path.exists(path):
				os.makedirs(path)

		# run Mosaic on `site` mosaics 
		Mosaic(resolution, '.', verbose=verbose, bkmodel=bkmodel)

class Field:
	"""
	Mosaic manager for `Montage`. 
	"""
	def __init__(self, center, sides, grid, subgrid, **kwargs):
		"""
		Initialize a Field centered on `center` (array-like of length two with
		right ascension and declination in decimal degrees) and side lengths 
		`sides` (array-like of length two in units of decimal degrees in 
		right ascension and declination respectively). `grid` (array-like of
		length two) specifieds the grid layout (e.g., (4,4) says 4x4) to 
		subdivide the Field. The `subgrid` is the same but pertaining to the 
		further subdivision of each SubField into a grid layout of `sites`.
		See SubField class.

		kwargs = {
				'verbose': True   , # display message, progress
				'survey' : 'DSS'  , # 2MASS, SDSS, DSS 
				'band'   : 'DSS2B'  # filter `band` pertaining to `survey`
			}
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs,
				{
					'verbose': True   , # display message, progress
					'survey' : 'DSS'  , # 2MASS, SDSS, DSS 
					'band'   : 'DSS2B'  # filter `band` pertaining to `survey`
				})

			# function parameter assignments 
			verbose = options('verbose')
			survey  = options('survey')
			band    = options('band')
		
		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in Field().')

		# available bands for each survey 
		bands = {
				# Two-Micron All-Sky Survey
				'2MASS': ['J', 'H', 'K'],
					
				# Sloan Digital Sky Survey
				'SDSS': ['U','G','R','I','Z'],

				# STScI Digitized Sky Survey 
				'DSS': ['DSS1B','DSS1R','DSS2B','DSS2R','DSS2IR','Quick-V']
			}

		# check for appropriate survey, band 
		if survey not in bands:
			raise MontageError('`{}` was not a recognized survey for '
			'the Montage Python module.'.format(survey))
		if band not in bands[survey]:
			raise MontageError('`{}` was not a recognized filter band '
			'for the `{}` survey.'.format(band, survey))

		# check arguments 
		if ( not hasattr(center, '__iter__') or 
			not hasattr(sides, '__iter__') or not hasattr(grid, '__iter__') or 
			not hasattr(subgrid, '__iter__')):
			raise MontageError('Field() expects array-like arguments.')
		if ( len(center) != 2 or len(sides) != 2 or len(grid) != 2 or 
				len(subgrid) != 2 ):
			raise MontageError('Field() expects arguments to be length two.')

		if verbose:
			stdout.write('\n Setting up {:d}x{:d} Field around ({:.2f}, '
				'{:.2f}) ... '.format(grid[0], grid[1], center[0], center[1]))
			stdout.flush()

		# SolveGrid()
		self.ra_centers, self.dec_centers = SolveGrid(sides, grid)
		# relative to Field `center`:
		self.ra_centers  += center[0]
		self.dec_centers += center[1] 

		# side lengths for SubField`s
		self.sub_sides = ( sides[0] / grid[0], sides[1] / grid[1] )

		# set current directory to `Field directory`.
		self.location = os.path.abspath('.')

		# name SubField directories 
		self.folders = [os.path.join(self.location,'SubField-{:03d}'.format(a+1))
			for a in range(len(self.ra_centers) * len(self.dec_centers)) ]

		# create SubField directories
		for folder in self.folders:
			if not os.path.exists(folder):
				os.makedirs(folder)

		# zip together all ra, dec pairs
		self.sub_centers = [ (ra, dec) 
			for ra in self.ra_centers for dec in self.dec_centers ]

		# initialize empty list of SubField`s
		self.subfields = []

		# initalize all SubField`s
		for a, folder, sub_center in zip(range(len(self.folders)), self.folders, 
			self.sub_centers):

			if verbose:
				stdout.write('\n Initializing {} ... '.format(
					os.path.basename(folder)))
				stdout.flush()

			# change directories
			os.chdir(folder)

			# create SubField 
			self.subfields.append( SubField(sub_center, self.sub_sides,
				subgrid, survey = survey, band = band) )

			if verbose: stdout.write('done')

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def ArchiveList(self, **kwargs):
		"""
		Run `ArchiveList()` on all SubFields. The only keyword option is 
		`verbose` (default: True).
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in '
			'Field.ArchiveList().')

		if verbose:
			stdout.write('\n Running ArchiveList() on all SubFields ... ')
			stdout.write('\n ' + '-' * 70 + '\n')
			stdout.flush()

		# Run ArchiveList() on all SubField`s
		for a, subfield in enumerate(self.subfields):

			if verbose:
				stdout.write('\n Running ArchiveList() on {} ... '
					.format(os.path.basename(self.folders[a])))
				stdout.flush()
			
			# run ArchiveList()
			subfield.ArchiveList(verbose = False)

			if verbose: stdout.write('done')

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def ArchiveExec(self, **kwargs):
		"""
		Run `ArchiveExec()` on all SubFields. The only keyword option is 
		`verbose` (default: True).
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in '
			'Field.ArchiveExec().')

		if verbose:
			stdout.write('\n Running ArchiveExec() on all SubFields ... ')
			stdout.write('\n' + '-' * 70 + '\n')
			stdout.flush()

		# Run ArchiveExec() on all SubField`s
		for a, subfield in enumerate(self.subfields):

			if verbose:
				stdout.write('\n Running ArchiveExec() on `{}`: {} / {} ... '
					.format(os.path.basename(self.folders[a]), a + 1, 
					len(self.folders)))
				stdout.flush()
			
			# run ArchiveList()
			subfield.ArchiveExec(verbose = False)

			if verbose: stdout.write('done')

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def Build(self, resolution, **kwargs):
		"""
		Build(resolution, **kwargs):

		Run the build process for all SubFields in this Field. See the 
		documentation for Montage.Mosaic() and SubField.Build().

		kwargs = {
				'verbose': True, # display messages, progress
				'bkmodel': True  # run background modelling procedures.
			}
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True, # display messages, progress 
					'bkmodel': True  # run background modelling procedure
				})

			# function parameter assignments 
			verbose = options('verbose')
			bkmodel = options('bkmodel')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in '
			'Field.ArchiveExec().')

		if verbose:
			stdout.write('\n Running Build() on all SubField`s ... ')
			stdout.write('\n' + '=' * 70 + '\n')
			stdout.write(       '=' * 70 + '\n')
			stdout.flush()

		# Run ArchiveExec() on all SubField`s
		for a, subfield in enumerate(self.subfields):

			if verbose:
				stdout.write('\n Running Build() on `{}`: {} / {} ... '
					.format(os.path.basename(self.folders[a]), a + 1, 
					len(self.folders)))
				stdout.write('\n' + '=' * 70 + '\n')
				stdout.flush()
			
			# run ArchiveList()
			subfield.Build(verbose = verbose, bkmodel = bkmodel)

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def Collect(self, **kwargs):
		"""
		Run Collect() on all SubFields of this Field.
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True, # display messages, progress 
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in '
			'Field.ArchiveExec().')

		if verbose:
			stdout.write('\n Running Collect() on all SubField`s ... ')
			stdout.write('\n' + '=' * 70 + '\n')
			stdout.write(       '=' * 70 + '\n')
			stdout.flush()

		# Run ArchiveExec() on all SubField`s
		for a, subfield in enumerate(self.subfields):

			if verbose:
				stdout.write('\n Running Collect() on `{}`: {} / {} ... '
					.format(os.path.basename(self.folders[a]), a + 1, 
					len(self.folders)))
				stdout.write('\n' + '=' * 70 + '\n')
				stdout.flush()
			
			# run ArchiveList()
			subfield.Collect( verbose=verbose )

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def Merge(self, resolution, **kwargs):
		"""
		Merge(resolution, **kwargs):

		Run Merge() on all SubFields of this Field. See SubField.Merge().
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True, # display messages, progress 
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in '
			'Field.ArchiveExec().')

		if verbose:
			stdout.write('\n Running Merge() on all SubField`s ... ')
			stdout.write('\n' + '=' * 70 + '\n')
			stdout.write(       '=' * 70 + '\n')
			stdout.flush()

		# Run ArchiveExec() on all SubField`s
		for a, subfield in enumerate(self.subfields):

			if verbose:
				stdout.write('\n Running Merge() on `{}`: {} / {} ... '
					.format(os.path.basename(self.folders[a]), a + 1, 
					len(self.folders)))
				stdout.write('\n' + '=' * 70 + '\n')
				stdout.flush()
			
			# run ArchiveList()
			subfield.Merge( verbose=verbose )

		if verbose:
			stdout.write('\n')
			stdout.flush()

	def Finalize(self, resolution, **kwargs):
		"""
		Finalize(resolution, **kwargs):

		Collect all SubField/master mosaics into a single folder and 
		run Mosaic() on them for a single final image.
		"""
		try:

			# function parameter defaults 
			options = Options( kwargs, 
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments 
			verbose = options('verbose')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment on Field.Finalize().')

		# relocate to Field directory 
		os.chdir(self.location)

		# create master directory 
		master_dir = os.path.join(self.location, 'master')
		if not os.path.exists(master_dir):
			os.makedirs(master_dir)

		# create master directory sub-structure
		for subdir in ['images','projected','differences','corrected','final']:
			path = os.path.join(master_dir,subdir)
			if not os.path.exists(path):
				os.makedirs(path)

		
		# collect master mosaics 
		for a, folder in enumerate(self.folders):

			if verbose:
				stdout.write('\n Copying image from {} ... '.format(
					os.path.basename(folder)))
				stdout.flush()
				
			# SubField mosaic 
			image = os.path.join(folder,'master/final/mosaic.fits')

			# ensure Merge() was run 
			if not os.path.exists(image):
				raise MontageError('Missing `master` mosaic image in `{}`.'
				.format(os.path.basename(folder)))

			# copy image to Field `master` image directory
			sh.copy(image, os.path.join(master_dir, 'images/mosaic-{}.fits'
				.format(a + 1)))

			if verbose: stdout.write('done')

		if verbose:
			stdout.write('\n')
			stdout.flush()

		# change directories to `master`
		os.chdir(master_dir)

		# run Mosaic() on all SubField `master`s
		Mosaic(resolution, '.', verbose=verbose, bkmodel=bkmodel)
