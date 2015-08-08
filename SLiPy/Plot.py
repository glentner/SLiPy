# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# AstroPython/Astro/Plot.py
"""
Plotting facility for astronomy.
"""
import matplotlib as mpl
from matplotlib import pyplot as plt

from .. import SlipyError
from ..Framework.Options import Options, OptionsError
from . import Fits, Spectrum

# mpl.rcParams['figure.facecolor'] = 'w'
plt.ion()

class PlotError(SlipyError):
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
        Assign `options` in `kwargs` and initialize the figure.
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
            self.ylimits  = []
            self.gridv    = None
            self.yargs    = []
            self.xargs    = []
            self.largs    = []
            self.targs    = []
            self.xkwargs  = {}
            self.ykwargs  = {}
            self.tkwargs  = {}
            self.lkwargs  = {}
            self.txargs   = []
            self.txkwargs = []

            if type(spectrum) is not Spectrum.Spectrum:
                raise PlotError('Splot expects type `Spectrum`!')

            # data in `list` allows for overplotting
            self.data = [ spectrum.data ]
            self.wave = [ spectrum.wave ]

            if spectrum.data.unit:
            	self.yargs = [str(spectrum.data.unit)]

            if spectrum.wave.unit:
            	self.xargs = [str(spectrum.wave.unit)]

            # `name` always retains `label`
            self.name   = self.options('label')

            # `label` and `marker`s same as data
            self.label  = [ self.options('label')  ]
            self.marker = [ self.options('marker') ]

            # set x limits to the data
            self.xlimits = [ spectrum.wave.min().value,
            	spectrum.wave.max().value ]

        except OptionsError as err:
            print(' --> OptionsError:', err.msg)
            raise PlotError('Failed to construct Splot.')

        # initialize the figure
        self.fig = plt.figure("Spectral-Plot (SLiPy)")
        self.ax  = self.fig.add_subplot(111)
        self.draw()

    def xlim(self, xmin, xmax ):
        """
        Handle to pyplot.xlim
        """
        self.xlimits = [ xmin, xmax ]
        self.ax.set_xlim(xmin, xmax)

    def ylim(self, ymin, ymax ):
        """
        Handle to pyplot.ylim
        """
        self.ylimits = [ ymin, ymax ]
        self.ax.set_ylim(ymin, ymax)

    def xlabel(self, *args, **kwargs ):
        """
        x axis label.
        """
        self.xargs = args
        self.xkwargs = kwargs
        self.ax.set_xlabel( *args, **kwargs )

    def ylabel(self, *args, **kwargs ):
        """
        y axis label.
        """
        self.yargs = args
        self.ykwargs = kwargs
        self.ax.set_ylabel( *args, **kwargs )

    def title(self, *args, **kwargs ):
        """
        title for plot.
        """
        self.targs = args
        self.tkwargs = kwargs
        self.ax.set_title( *args, **kwargs )

    def legend(self, *args, **kwargs):
        """
        legend for plot.
        """
        self.largs = args
        self.lkwargs = kwargs
        plt.legend( *args, **kwargs )

    def text(self, *args, **kwargs):
        """
        display text over plot.
        """
        self.txargs.append( args )
        self.txkwargs.append( kwargs )
        self.ax.text( *args, **kwargs )

    def txtclear(self):
        """
        Clear all `text` from figure.
        """
        self.txargs = []
        self.txkwargs = []
        self.draw()

    def markers(self, *args):
        """
        Reassign the values for the `marker`s in the figure. The number
        of arguments must equal the number of spectra in the figure. This
        starts out as one, but will increase for ever SPlot.overlay().
        """
        if len(args) != len(self.data):
            raise PlotError('{} arguments were given but there are {} '
            'spectra plotted in this figure!'.format(len(args), len(self.data)))

        for a, mark in enumerate(args):
            if type(mark) is not str:
                raise PlotError('Arguments given to SPlot.markers() must be '
                '{} but argument #{} was {}'.format(type(''), a+1, type(mark)))

        self.marker = list(args)

    def __build(self, picker = False):
        """
        Make the plot.
        """

        if picker:
            self.restore()
            self.ax.plot(self.wave[0], self.data[0], self.marker[0],
                label = self.label[0], picker = True)

        else:
            for x, y, m, l in zip(self.wave, self.data, self.marker, self.label):
            	self.ax.plot(x, y, m, label=l)

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
        	    self.ax.text( *args, **kwargs )

        if self.xlimits:
            self.xlim( *self.xlimits )

        if self.ylimits:
            self.ylim( *self.ylimits )

        if self.gridv:
            self.grid(self.gridv)

        if self.usetex:
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')

    def refresh(self):
        """
        pyplot.draw()
        """
        plt.draw()

    def draw(self, picker = False):
    	"""
    	Re-build the plot
    	"""
    	self.ax.clear()
    	self.__build(picker = picker)
    	plt.draw()

    def close(self):
        """
        Close the plot.
        """
        plt.close("Spectral-Plot (SLiPy)")

    def grid(self, value):
        """
        Show grid on plot.
        """
        self.gridv = value
        plt.grid(value)

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

    def yoffset(self, value):
        """
        Toggle the offset for the y axis
        """
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(value)

    def tight_layout(self):
        """
        pyplot.tight_layout()
        """
        plt.tight_layout()

    def overlay(self, *splots ):
        """
        Overlay (add) spectra to this plot from other `splots`.
        """

        for a, plot in enumerate(splots):

            # check data type
            if type(plot) is not SPlot:
                raise PlotError('Splot.overlay expects '
                'type Splot! (from argument {})'.format(a))

            # add data
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
