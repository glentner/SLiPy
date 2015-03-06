# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/Montage.py 
"""
This module makes use of the `Montage` mosaic tools from caltech, see:
http://montage.ipac.caltech.edu/

The user should have Montage`s executables available on their path.
"""

import os, subprocess, shutil as sh, numpy as np
from numbers import Number
from ..Framework.Options import Options, OptionsError
from ..Framework.Display import Monitor, DisplayError


class MontageError(Exception):
	"""
	Exception specific to the Montage module
	"""
	pass

class SubField:
	"""
	SubField( center, region=(6,6), scale=3, **kwargs ):
	"""
	def __init__(self, center, region = (6,6), scale = 3, **kwargs ):
		"""
		Solve for the the centers if multiple `sites` given a `center` 
		location (tuple of length two containing the right ascension and
		declination in decimal degrees), the size of the `region` (tuple
		of length two containing the horizontal and vertical components
		in degrees), and the `scale` of the grid space (resolution).

		The object builds the commands needed for Montage using the 
		`mArchiveList` command on the `sites`. `survey` and `band` must be
		specified in keyword arguments (see `Montage`s documentation).
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
			if not hasattr( center, '__iter__'):
				raise MontageError('Montage.SubField() expects a tuple of '
				'length two for `center` argument.')
			if not hasattr( region, '__iter__'):
				raise MontageError('Montage.SubField() expects a tuple of '
				'length two for `region` argument.')
			if len(center) != 2:
				raise MontageError('Montage.SubField() expects `center` '
				'argument to have exactly two elements.')
			if len(region) != 2:
				raise MontageError('Montage.SubField() expects `region` '
				'argument to have exactly two elements.')
			if not isinstance( scale, Number):
				raise MontageError('Montage.SubField() expects `scale` '
				'argument to be a number.')

			# grid `site` centers in the horizontal axis
			ra_left_site    = -region[0]/2 + scale / 2
			ra_right_site   =  region[0]/2 - scale / 2 
			ra_site_centers = np.linspace( ra_left_site, ra_right_site,
				(ra_right_site - ra_left_site) / scale + 1 )

			# grid `site` centers in the vertical axis 
			dec_bottom_site  = -region[1]/2 + scale / 2
			dec_top_site     =  region[1]/2 - scale / 2
			dec_site_centers = np.linspace( dec_bottom_site, dec_top_site,
				(dec_top_site - dec_bottom_site) / scale + 1)
			
			# relative to SubField `center`:
			ra_site_centers  += center[0]
			dec_site_centers += center[1] 

			# record number of `site`s along each axis
			self.num_ra_sites  = len(ra_site_centers)
			self.num_dec_sites = len(dec_site_centers) 

			# build arguments for subprocess call
			self.ArchiveList = [
					['mArchiveList', survey, band, 
						'{:.2f} {:.2f}'.format(ra_site, dec_site), str(scale), 
						str(scale), 'remote.tbl']
					for ra_site in ra_site_centers for dec_site in dec_site_centers
				]
			
		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assignment in Subfield().')

	def archive(self, **kwargs):
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

			# create directory trees
			subdirs = [ 
					'SubField_{}{}/raw'.format(a, b) 
					for a in range(self.num_ra_sites) 
					for b in range(self.num_dec_sites) 
				]
			
			for tree in subdirs:
				os.makedirs(tree)

			if verbose:
				display = Monitor()
				nargs   = len(self.ArchiveList)
				print(' Running `mArchiveList` at {} sites ... '.format(nargs))

			for a, command in enumerate(self.ArchiveList):
				# navigate to `raw` directory 
				os.chdir( subdirs[a] )
				# submit subprocess call
				output = subprocess.check_output(command).decode('utf-8')
				# check output for success
				if 'ERROR' in output or 'count="0"' in output:
					raise MontageError('Failed archive list from archive() '
					'(command: {}), (output: {}).'.format(command, output))
				# display progress 
				if verbose: display.progress(a, nargs)
				# move up directory tree 
				os.chdir('../../')

			# erase progress bar
			if verbose: display.complete()

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise MontageError('Failed keyword assigned from SubField.exec().')

		except DisplayError as err:
			print(' --> DisplayError:', err)
			raise MontageError('Display.Monitor() failure in SubField.exec().')
