# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/Montage.py 
"""
This module makes use of the `Montage` mosaic tools from caltech, see:
http://montage.ipac.caltech.edu/

The user should have Montage`s executables available on their path.
"""

import os, shutil as sh, numpy as np 
from subprocess import check_output as call
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

class SubField:
	"""
	SubField( center, region=(6,6), scale=3, **kwargs ):
	"""
	def __init__(self, center, sides, grid, **kwargs ):
		"""
		Create `site` grid for SubField.
		"""
		try:
			# function parameter defaults
			options = Options( kwargs,
				{
					'survey': '', # DSS, SDSS, 2MASS
					'band'  : '', # filter for `survey`, see `bands` dictionary
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
			if not survey or not band:
				raise MontageError('Both the `survey` and `band` must be '
				'specified in keyword arguments for SubField().')
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

			# new tree structure
			self.folders = [ os.path.join( os.path.abspath('.'),
					'SubField_{}{}'.format(a + 1, b + 1) )
					for a in range(self.num_ra_sites) 
					for b in range(self.num_dec_sites) 
				]
		
			# initialize folder structure with `images` directory
			for folder in self.folders:
				abspath = os.path.join(folder,'images')
				if not os.path.exists(abspath):
					os.makedirs(abspath)

			if verbose:
				display = Monitor()
				nargs   = len(self.archive_command_list)
				print(' Running `mArchiveList` at {} sites ... '.format(nargs))

			for a, command in enumerate(self.archive_command_list):
				# navigate to `raw` directory 
				os.chdir( os.path.join(self.folders[a],'images') )
				# submit subprocess call
				output = call(command).decode('utf-8')
				# check output for success
				if 'ERROR' in output or 'count="0"' in output:
					raise MontageError('Failed archive list from archive() '
					'(command: {}), (output: {}).'.format(command, output))
				# display progress 
				if verbose: display.progress(a, nargs )

			# erase progress bar
			if verbose: display.complete()

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigned from SubField.exec().')

		except DisplayError as err:
			print(' --> DisplayError:', err)
			raise MontageError('Display.Monitor() failure in SubField.'
			'ArchiveList().')

	def ArchiveExec(self, **kwargs):
		"""
		Run `mArchiveExec` on each site in the SubField.
		"""
		try:
			# function parameter defaults 
			options = Options( kwargs,
				{
					'verbose': True # display messages, progress
				})

			# function parameter assignments
			verbose = options('verbose')

			# initialize display 
			if verbose: 
				display  = Monitor()
				nfolders = len(self.folders)
				print(' Running `mArchiveExec` on {} sites ... '
					.format(nfolders))

			for a, folder in enumerate(self.folders):
				# navigate to site folder 
				os.chdir( os.path.join(folder,'images') )
				# run `mArchiveExec`
				output = call(['mArchiveExec','remote.tbl']).decode('utf-8')
				# check for errors
				if 'ERROR' in output:
					raise MontageError('Failed `mArchiveExec` in folder `{}`.'
					'Output: {}'.format(folder, output))
				# display progress 
				if verbose: display.progress(a, nfolders)

			# erase progress bar 
			if verbose: display.complete()

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigment in ArchiveExec().')

		except DisplayError as err:
			print(' --> DisplayError:', err)
			raise MontageError('Display.Monitor failure in SubField.'
			'ArchiveExec().')

	def Build(self, res, **kwargs):
		"""
		Run the build process for the `sites` in this SubField. `res`olution 
		argument is the number of pixels per degree desired.
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
	
		# setup folder structure 
		if verbose: 
			stdout.write(' Setting up folder structure ... ')
			stdout.flush()
		
		for subdir in ['corrected','projected','differences','final']:
			for folder in self.folders:
				abspath = os.path.join(folder,subdir)
				if not os.path.exists(abspath):
					os.makedirs(abspath)
		
		if verbose: 
			stdout.write('done\n')
			stdout.flush()

		# build mosaics at all sites 
		for a, folder in enumerate(self.folders):
			
			# change directories 
			os.chdir(folder)

			# display message
			if verbose: 
				stdout.write(' Building mosaic for SubField-site: {} / {} ... \n'
				.format(a + 1, len(self.folders)))
				stdout.write( '-' * 70 + '\n' )
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
			output = call(['mMakeHdr','-p','{}'.format(1 / res),
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


