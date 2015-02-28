#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/Astro/Plot.py 
"""
Plotting facility for astronomy. 
"""
import matplotlib as mpl
from matplotlib import pyplot as plt 
from Python import BaseError
from Python.General.Options import *
from Python.Astro import Fits

mpl.rcParams['figure.facecolor'] = 'w'
plt.ion()

class PlotError(BaseError):
	"""
	Exception specific to Plot module.
	"""
	pass

class SPlot:
	"""
	SPlot( spectra, **kwargs )

	Spectrum Plot - Plot the data in `spectra`.
	"""
	def __init__(self, spectra, **kwargs):
		"""
		Assign `options` in `kwargs` and initialize the plot.
		"""
		try:
			# available options
			self.options = Options( kwargs,
				{
					'marker': 'b-'        , # marker for plot 
					'label' : ''          , # label for data 
					'xlabel': 'Wavelength', # label for x-axis
					'ylabel': ''          , # label for y-axis
					'title' : ''          , # title of plot 
					'labelpad': 10        , # label pad for x,y axis 
					'fontsize': 12        , # font size
					'usetex': False         # pdflatex setting
				})

			# assign options 
			self.xlabel   = self.options('xlabel')
			self.ylabel   = self.options('ylabel')
			self.title    = self.options('title')
			self.labelpad = self.options('labelpad')
			self.fontsize = self.options('fontsize')
			self.usetex   = self.options('usetex')

			if type(spectra) is not Fits.Spectra:
				raise PlotError('Splot expects type Fits.Spectra!')

			# data in `list` allows for overplotting
			self.data = [ spectra.data ]
			
			if spectra.wavecal:
				self.wavecal = True
				self.wave    = [ spectra.wave ]

			else: self.wave = [ 
				np.arange( np.shape(spectra.data)[0] ) ]

			self.label  = [ self.options('label')  ]
			self.marker = [ self.options('marker') ]

			# set x limits to the data 
			if self.wavecal: self.xlimits = [ 
				spectra.wave.min(), spectra.wave.max() ]

			else: self.xlimits = [
					0, np.shape(spectra.data)[0] ]

		except OptionsError as err:
			print(' --> OptionsError:', err.msg)
			raise PlotError('Failed to construct Splot.')

	def xlim(self, xmin, xmax ):
		"""
		Handle to pyplot.xlim
		"""
		plt.xlim(xmin, xmax)

	def ylim(self, ymin, ymax ):
		"""
		Handle to pyplot.ylim
		"""
		plt.ylim(ymin, ymax)

	def __build(self):
		"""
		Make the plot.
		"""
		for x, y, m, l in zip(
			self.wave, self.data, self.marker, self.label):
			
			plt.plot(x, y, m, label=l)
				
		plt.xlabel(self.xlabel, labelpad = self.labelpad,
			fontsize = self.fontsize )

		plt.ylabel(self.ylabel, labelpad = self.labelpad,
			fontsize = self.fontsize )
				
		plt.title( self.title, fontsize = self.fontsize )

		self.xlim( *self.xlimits )

		if self.usetex:
			plt.rc('text', usetex=True)
			plt.rc('font', family='serif')
		
	def draw(self):
		"""
		Remake plot and run pyplot.draw 
		"""
		plt.clf()
		self.__build()
		plt.draw()

	def show(self):
		"""
		Show the plot 
		"""
		plt.show()

	def clf(self):
		"""
		Clear the plot.
		"""
		plt.clf()

	def close(self):
		"""
		Close the plot.
		"""
		plt.close()

	def gridon(self):
		"""
		Show grid on plot.
		"""
		plt.grid(True)
		plt.draw()

	def gridoff(self):
		"""
		Turn off grid for plot.
		"""
		plt.grid(False)
		plt.draw()

	def legend(self, *args, **kwargs):
		"""
		legend(*args, **kwargs):

		Same as pyplot.legend.
		"""
		plt.legend(*args, **kwargs)
	
	def save(self, filename):
		"""
		Save plot to `filename`. Must have extension for formatting.
		"""
		if type(filename) is not str:
			raise PlotError('`filename` should be of type str.')
		if len(filename.split('.')) < 2:
			raise PlotError('`filename` needs an extension.')

		plt.savefig(filename, format=filename.split('.')[-1])

	def overlay(self, *splots ):
		"""
		Overlay (add) spectra to this plot from other `splots`.
		"""
		# check data type 
		for a, plot in enumerate(splots):
			if type(plot) is not Splot:
				raise PlotError('Splot.overlay expects '
					'type Splot! (from argument {})'.format(a))
			
			if self.wavecal and not plot.wavecal:
				raise PlotError('Original spectrum is wavelength '
					'calibrated but Splot {} in arguments was not!'
					.format(a))

			self.data   += plot.data 
			self.wave   += plot.wave 
			self.marker += plot.marker 
			self.label  += plot.label 

	def restore(self):
		"""
		Restore self.data and self.wave from possible `overlay`s. 
		"""
		self.data   = [ self.data[0]   ]
		self.wave   = [ self.wave[0]   ]
		self.marker = [ self.marker[0] ]
		self.label  = [ self.label[0]  ]


			
