#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/Astro/Plot.py 
"""
Plotting facility for astronomy. 
"""
import matplotlib as mpl
from matplotlib import pyplot as plt 

from ..Framework.Options import Options, OptionsError
from . import Fits, DataType

mpl.rcParams['figure.facecolor'] = 'w'
plt.ion()

class PlotError(Exception):
	"""
	Exception specific to Plot module.
	"""
	pass

class SPlot:
	"""
	SPlot( spectrum, **kwargs )

	Spectrum Plot - Plot the data in `spectrum`.
	"""
	def __init__(self, spectrum, **kwargs):
		"""
		Assign `options` in `kwargs` and initialize the plot.
		"""
		try:
			# available options
			self.options = Options( kwargs,
				{
					'marker': 'b-'          , # marker for plot 
					'label' : 'unspecified' , # label for data 
					'usetex': False           # pdflatex setting
				})

			# assign options 
			self.usetex   = self.options('usetex')
			self.ylimits  = None
			self.gridv    = None
			self.yargs    = None
			self.xargs    = None
			self.largs    = None
			self.targs    = None
			self.xkwargs  = None
			self.ykwargs  = None 
			self.tkwargs  = None
			self.lkwargs  = None
			self.txargs   = []
			self.txkwargs = []

			if type(spectrum) is not DataType.Spectrum:
				raise PlotError('Splot expects type Fits.Spectra!')

			# data in `list` allows for overplotting
			self.data = [ spectrum.data ]
			
			if spectrum.wavecal:
				self.wavecal = True
				self.wave    = [ spectrum.wave ]

			else: self.wave = [ 
				np.arange( np.shape(spectrum.data)[0] ) ]

			# `name` always retains `label`
			self.name   = self.options('label')
			# `label` and `marker`s same as data
			self.label  = [ self.options('label')  ]
			self.marker = [ self.options('marker') ]

			# set x limits to the data 
			if self.wavecal: self.xlimits = [ 
				spectrum.wave.min(), spectrum.wave.max() ]

			else: self.xlimits = [ 0, len(spectrum.data) ]

		except OptionsError as err:
			print(' --> OptionsError:', err.msg)
			raise PlotError('Failed to construct Splot.')

	def xlim(self, xmin, xmax ):
		"""
		Handle to pyplot.xlim
		"""
		self.xlimits = [ xmin, xmax ]
		plt.xlim(xmin, xmax)
		plt.draw()

	def ylim(self, ymin, ymax ):
		"""
		Handle to pyplot.ylim
		"""
		self.ylimits = [ ymin, ymax ]
		plt.ylim(ymin, ymax)
		plt.draw()

	def xlabel(self, *args, **kwargs ):
		"""
		x axis label.
		"""
		self.xargs = args 
		self.xkwargs = kwargs 
		plt.xlabel( *args, **kwargs )
		plt.draw()
	
	def ylabel(self, *args, **kwargs ):
		"""
		y axis label.
		"""
		self.yargs = args
		self.ykwargs = kwargs 
		plt.ylabel( *args, **kwargs )
		plt.draw()

	def title(self, *args, **kwargs ):
		"""
		title for plot.
		"""
		self.targs = args 
		self.tkwargs = kwargs
		plt.title( *args, **kwargs )
		plt.draw()

	def legend(self, *args, **kwargs):
		"""
		legend for plot.
		"""
		self.largs = args 
		self.lkwargs = kwargs 
		plt.legend( *args, **kwargs )
		plt.draw()

	def text(self, *args, **kwargs):
		"""
		display text over plot.
		"""
		self.txargs.append( args )
		self.txkwargs.append( kwargs )
		plt.text( *args, **kwargs )
		plt.draw()

	def txtclear(self):
		"""
		Clear all `text` from figure.
		"""
		self.txargs = []
		self.txkwargs = []
		self.draw()

	def __build(self):
		"""
		Make the plot.
		"""
		for x, y, m, l in zip(
			self.wave, self.data, self.marker, self.label):
			plt.plot(x, y, m, label=l)
				
		if self.xargs or self.xkwargs: 
			self.xlabel( *self.xargs, **self.xkwargs ) 
		
		if self.yargs or self.ykwargs: 
			self.ylabel( *self.yargs, **self.ykwargs )
		
		if self.targs or self.tkwargs: 
			self.title( *self.targs, **self.tkwargs )
		
		if self.largs or self.lkwargs: 
			self.legend( *self.largs, **self.lkwargs )
		
		if self.txargs or self.txkwargs:
			for args, kwargs in zip(self.txargs, self.txkwargs):
				plt.text( *args, **kwargs )

		if self.xlimits: 
			self.xlim( *self.xlimits )
		
		if self.ylimits: 
			self.ylim( *self.ylimits )
		
		if self.gridv: 
			self.grid(self.gridv)
		
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
		plt.draw()

	def close(self):
		"""
		Close the plot.
		"""
		plt.close()

	def grid(self, value):
		"""
		Show grid on plot.
		"""
		self.gridv = value
		plt.grid(value)
		plt.draw()

	def save(self, filename):
		"""
		Save plot to `filename`. Must have extension for formatting.
		"""
		if type(filename) is not str:
			raise PlotError('`filename` should be of type str.')
		if len(filename.split('.')) < 2:
			raise PlotError('`filename` needs an extension.')

		plt.savefig(filename, format=filename.split('.')[-1])

	def xoffset(self, value):
		"""
		Toggle the offset for the x axis
		"""
		plt.gca().get_xaxis().get_major_formatter().set_useOffset(value)
		plt.draw()
	
	def yoffset(self, value):
		"""
		Toggle the offset for the y axis 
		"""
		plt.gca().get_yaxis().get_major_formatter().set_useOffset(value)
		plt.draw()

	def overlay(self, *splots ):
		"""
		Overlay (add) spectra to this plot from other `splots`.
		"""
		# check data type 
		for a, plot in enumerate(splots):
			if type(plot) is not SPlot:
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

def desired( plot ):
	"""
	desired( plot ):

	Helper function for Iterate. Prompts user to keep `plot`;
	returns True or False.
	"""
	# draw the plot
	plot.draw()
	# prompt the user for input
	prompt = input('\r\033[K keep -> `{}` (y/[n]/x)?: '
			.format(plot.name)).strip()
	# insure valid response
	while True:
		if prompt not in ['y','n','','x']:
			# invalid input, prompt again
			print('\r\033[K `{}` was not a recognized response.'.format(prompt))
			prompt = input('\033[2A\r\033[K keep -> `{}` (y/[n]/x)?: '
					.format(plot.name)).strip()
		else: 
			# clear the error message
			print('\r\033[K\033[1A')
			break
	
	if prompt in ['n', '']:
		return False

	elif prompt in ['y']:
		return True

	else: 
		# the user input `x`
		raise KeyboardInterrupt('\r\033[K User exitted early, saving results.')
	

def Iterate( *plots, **kwargs ):
	"""
	Iterate( *plots, **kwargs ):

	Iterate thru `plots` to inspect data, the user marks `plots` of 
	interest. The function returns a list of `names` marked.
	"""
	try:
		options = Options( kwargs,
			{
				'keep' : 'name' # alternatively, `plot`
			})

		keep = options('keep')

		if keep not in ['name', 'plot']:
			raise PlotError('Iterate expects either `name` or `plot` for '
				'keyword argument `keep`.')

		# check input arguments
		for plot in plots:
			if not hasattr(plot, 'draw'):
				raise PlotError('Iterate expects objects to '
					'have a `draw` method.')
			if not hasattr(plot, 'name'):
				raise PlotError('Iterate expects objects to '
					'have a `name` method.')

		# clear some space
		print('\n')

		keepers = []
		for a, plot in enumerate(plots):
			print('\033[2A\r\033[K Showing plot {} of {} ... '
					.format(a, len(plots)) )
			if desired( plot ):
				if keep == 'name':
					keepers.append( plot.name )
				elif keep == 'plot':
					keepers.append( plot )

		return keepers

	except OptionsError as err:
		print(' --> OptionsError:', err.msg)
		raise PlotError('Failed to initialize Iterate.')

	except KeyboardInterrupt as x:
		print(x)
		return keepers
