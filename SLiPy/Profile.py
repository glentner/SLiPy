# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/SLiPy/Profile.py

"""
Profile fitting tasks for spectra.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate.interpolate import interp1d
from scipy.special import wofz as w

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import widgets

from astropy import units as u

from .. import SlipyError
from .Observatory import Observatory
from .Spectrum import Spectrum, SpectrumError
from .Plot import SPlot, PlotError
from ..Framework.Options import Options, OptionsError

from ..Algorithms.Functions import Gaussian, Lorentzian, Voigt, InvertedLorentzian
from ..Algorithms.KernelFit import KernelFit1D

class ProfileError(SlipyError):
	"""
	Exception specific to Profile module.
	"""
	pass

def Pick(event):
    """
    Used to hand selection events.
    """

    if isinstance(event.artist, Line2D):

        # get the x, y data from the pick event
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind

        # update the selection dictionary
        global selected
        selected['wave'].append(np.take(xdata,ind)[0])
        selected['data'].append(np.take(ydata,ind)[0])

        # display points as a visual aid
        plt.scatter(np.take(xdata,ind)[0], np.take(ydata,ind)[0],
            marker = 'o', s=75, c='r')
        plt.draw()

# empty selection dictionary is filled with the Select() function
selected = {'wave': [], 'data': []}

def Select(splot):
    """
    Select points from the `splot`. This should be of type SPlot
    (or it can optionally be a Spectrum type, for which a SPlot will be
    created). The splot will be rendered and the user clicks on the
    figure. When finished, return to the terminal prompt. A dictionary is
    returned with two entries, `wave` and `data`, representing the x-y
    locations selected by the user. This can always be retrieved later by
    accessing the module member `Profile.selected`.
    """
    if type(splot) is Spectrum:
        splot = SPlot(splot)

    elif type(splot) is not SPlot:
        raise ProfileError('Select() requires either a Spectrum or SPlot '
        'object as an argument')

    # reset the selection dictionary
    global selected
    selected = { 'wave': [], 'data': []}

    splot.draw(picker = True)

    splot.fig.canvas.mpl_connect('pick_event', Pick)

    input(' Press <Return> after making your selections ... ')
    return selected

def AutoFit(splot, function = InvertedLorentzian, params = None):
    """
    Given `splot` of type SPlot, the user selects two points on the
    spectrum and a parameterized function is fit (an inverted Lorentzian by
    default). Optionally, `splot` can be of type spectrum and a basic SPlot
    will be created for you. If the user gives an alternative `function`,
    `params` (parameters) must be provided. `params` is to be the first guess,
    `p0` given to scipy...curve_fit; the user can provide them expicitely,
    or in the form of functions with the template `function(xarray, yarray)`
    where `xarray` and `yarray` are the `wave` and `data` arrays extracted
    between the two points selected by the user.
    """

    print(' Please select four points identifying the spectral line.')
    print(' Outer points mark the domain of the line.')
    print(' Inner points mark the sample of the line to fit.')

    # make selections
    selected = Select(splot)

    if len( selected['wave'] ) != 4:
        raise ProfileError('Exactly 4 locations should be selected for '
        'the Profile.Fit() routine!')

    # order the selected wavelength locations
    points = selected['wave']
    points.sort()

    # get data and wavelength arrays
    if type(splot) is SPlot:
        wave, data = splot.wave[0], splot.data[0]
    else:
        wave, data = splot.wave, splot.data

    # extract the domains `selected` (i.e., the sample interval)
    x_inner = wave[np.where(np.logical_and(points[1] < wave, wave < points[2]))]
    x_outer = wave[np.where(np.logical_and(points[0] < wave, wave < points[3]))]
    y_inner = data[np.where(np.logical_and(points[1] < wave, wave < points[2]))]
    y_outer = data[np.where(np.logical_and(points[0] < wave, wave < points[3]))]

    # y_inner = data[ wave[ wave < points[2] ] > points[1] ]
    # x_inner = wave[ wave[ wave < points[2] ] > points[1] ]

    if function.__name__ == 'InvertedLorentzian':
        # First guess for default behavior
        params = [ y_inner.min().value, x_inner[ y_inner.argmin() ].value,
            (x_inner[-1] - x_inner[0]).value / 2]

    elif not params:
        # the user gave a different function but not parameters!
        raise ProfileError('The user must provide `params` when giving an '
        'alternative `function` in Profile.Fit()!')

    else:

        if not hasattr(params, '__iter__'):
            raise ProfileError('`params` must be an iterable type in '
            'Profile.Fit()!')

        try:
            for a, parameter in enumerate(params):
                if type(parameter) is type(InvertedLorentzian):
                    # replace parameter with function evaluation
                    params[a] = parameter(x_inner, y_inner)

        except TypeError as err:
            print(' --> TypeError:', err)
            raise ProfileError('Profile.Fit() failed to call user functions '
            'correctly in `params`!')

    # fit a parameterized curve
    coeff, var_matrix = curve_fit(
            function,      # function to call, default is InvertedLorentzian
            x_inner.value, # domain (without units)
            y_inner.value, # data (without units)
            p0 = params    # list of parameters to try as a first guess
        )

    # display visual aids ...
    # evaluation of the fit profile over larger domain
    plt.plot(x_outer, function(x_outer.value, *coeff) * y_outer.unit,
        'b--', linewidth = 4)
    plt.plot(x_inner, function(x_inner.value, *coeff) * y_inner.unit,
        'r-', linewidth = 4)

    # return the larger domain, evaluated and of type Spectrum
    return Spectrum(function(x_outer.value, *coeff) * y_outer.unit, x_outer)


def Extract(splot, kernel = Gaussian, **kwargs):
    """
    Select locations in the `splot` figure, expected to be of type SPlot.
    Exactly four points should be selected. These are used to extract a
    line profile from the spectrum plotted in the splot figure. The inner
    section is used for the line, and the outer selection is used to model
    the continuum; these, respectively, and both returned as Spectrum objects.
    The gap is jumped using 1D interpolation (scipy...interp1d).
    """
    try:

        options = Options( kwargs, {
            'kind'      : 'cubic' , # given to scipy...interp1d for continuum
            'bandwidth' : 0.1*u.nm, # user should provide this!
            'rms'       : False     # measure the RMS of the line, continuum
        })

        kind      = options('kind')
        bandwidth = options('bandwidth')
        rms       = options('rms')

    except OptionsError as err:
        print(' --> OptionsError:', err)
        raise ProfileError('Unrecognized option given to Extract()!')

    print(' Please select four points identifying the spectral line.')
    print(' Outer intervals sample the continuum.')
    print(' Center interval contains the line.')

    # make selections
    selected = Select(splot)

    if len( selected['wave'] ) != 4:
        raise ProfileError('Exactly 4 locations should be selected for '
        'the profile modeling to work!')

    # order the selected wavelength locations
    wave = selected['wave']
    wave.sort()

    # create `line` profile
    xl = splot.wave[0].copy()
    yl = splot.data[0].copy()
    yl = yl[ xl[ xl < wave[2] ] > wave[1] ]
    xl = xl[ xl[ xl < wave[2] ] > wave[1] ]
    line = Spectrum(yl, xl)

    # extract continuum arrays
    xc = splot.wave[0].copy()
    yc = splot.data[0].copy()
    # inside outer-most selections
    yc = yc[ xc[ xc < wave[3] ] > wave[0] ]
    xc = xc[ xc[ xc < wave[3] ] > wave[0] ]
    # keep wavelengths whole domain for later
    xx = xc.copy()
    yy = yc.copy()
    # but not the line
    yc = yc[np.where(np.logical_or(xc < wave[1], xc > wave[2]))]
    xc = xc[np.where(np.logical_or(xc < wave[1], xc > wave[2]))]

    # use `kernel smoothing` to model the continuum
    model = KernelFit1D(xc, yc, kernel = kernel, bandwidth = bandwidth)

    # interpolate to cross `the gap`
    interp = interp1d(xc, model.mean(xc), kind = kind)

    # continuum model outside the line
    cont_outside = interp(xc)

    # continuum inside the line
    cont_inside = interp(xl)

    # continuum model for whole domain
    cont_domain = interp(xx)

    # build a spectrum from the arrays
    continuum = Spectrum(cont_domain, xx)

    # display visual aid
    plt.plot(xx, cont_domain, 'r--', linewidth = 3)
    plt.fill_between(xl, yl, cont_inside, color = 'blue', alpha = 0.25)
    plt.draw()

    if not rms:
        return line, continuum

    cont_rms = np.sqrt( np.sum( (cont_outside * yc.unit - yc)**2 ) / len(yc) )
    line_rms = cont_rms * np.sqrt(cont_inside / yl.value)

    return line, continuum, line_rms

class FittingGUI:
	"""
	Graphical Interface (keeps references alive) for fitting analytical
	profiles to spectral data.
	"""
	def __init__(self, fig, function='Voigt', **kwargs):
		"""
		Build widget elements based on a spectrum plotted in `fig` (type SPlot).
		`fig` may also be a Spectrum object (for which a SPlot will be
		created).

		Make initial guesses at parameters, build widgets, render the figure.
		"""
		try:

			options = Options( kwargs, {

					'kind'      : 'cubic' , # given to interp1d for continuum
					'bandwidth' : 0.1*u.nm, # user should provide this!
				})

			kind      = options('kind')
			bandwidth = options('bandwidth')

		except OptionsError as err:
			print(' --> OptionsError:', err)
			raise ProfileError('Unrecognized option given to Extract()!')

		if function not in ['Voigt', 'Lorentzian', 'Gaussian']:
			raise ProfileError('The only currently implemented functions '
			'for fitting profiles are `Voigt`, `Lorentzian` and `Gaussian`!')

		# grab and/or create the SPlot
		if isinstance(fig, Spectrum):
			fig = SPlot(fig, marker='k-', label='spectrum')

		elif not isinstance(fig, SPlot):
			raise ProfileError('FittingGUI() expects the `fig` argument to '
			'either be a Spectrum or a SPlot!')

		print('\n We need to extract the lines from the continuum before we '
		'begin the fitting process.')

		# extract the line, continuum, rms from the spectrum
		self.line, self.continuum, self.rms = Extract(fig, bandwidth=bandwidth,
			kind=kind, rms=True)

		print('\n Now select the peaks of each line to be fit.')
		print(' Initial guesses will be made for each line markered.')
		input(' Press <Return> after making your selections ... ')

		# grab all the selected points
		global selected
		points = np.array([
        	[ entry.value for entry in selected['wave'] ],
        	[ entry.value for entry in selected['data'] ]
			])

		# point pairs in ascending order by wavelength
		points = points[:, points[0].argsort()]

		# domain size of the line and the number of components
		self.domainsize = self.line.wave.value[-1] - self.line.wave.value[0]
		self.numlines   = np.shape(points)[1] - 4

		if self.numlines < 1:
			raise ProfileError('FittingGUI() expects at least one line to '
			'be selected for fitting!')

		# initial guesses for parameters given the line locations,
		# containing `L1`, `L2`, etc ... the values of which are themselves
		# dictionaries of parameters (e.g., `FWHM`, `Depth`, etc...) whose values
		# are the initial guesses for those parameters given the location of
		# it`s peak and the number of peaks within the `line`s domain.
		self.Params = {

				'L' + str(line + 1) : self.Parameterize(function, loc)
				for line, loc in enumerate(points[:, 2:-2].T)

			}

		# final spectral line profile is convolution of the LSP and the
		# line profile function. Gaussian on Gaussian, Gaussian on Lorentzian,
		# etc ... Set the functional form given requested profile function
		self.Evaluate = self.SetFunction(function)

		# grab the actual Figure object and it's axis to keep references
		self.fig = fig.fig
		self.ax  = fig.ax

		# refresh image, but keep any changes in axis limits
		fig.xlim( *fig.ax.get_xlim() )
		fig.ylim( *fig.ax.get_ylim() )
		fig.draw()

		# bring up plot to make room for sliders
		plt.subplots_adjust(bottom = 0.30)

		# resample continuum onto line domain
		self.continuum = self.continuum.copy()
		self.continuum.resample(self.line)

		# common domain for plotting, strip units, wavelengths from continuum
		self.x         = self.continuum.wave.value
		self.continuum = self.continuum.data.value

		# copy of continuum allows for adjustments
		self.y = self.continuum.copy()

		# save the mid-level of the continuum to avoid over-calculation
		self.continuum_level = self.continuum.mean()

		# add plots of each line, keep dictionary of handles
		self.Component = {

			line : plt.plot(self.x,
				self.y - self.Evaluate(self.x, **parameters), 'k--')[0]

			for line, parameters in self.Params.items()
		}

		# add plot of superposition of each component
		self.Combination, = plt.plot(self.x, self.y - self.SuperPosition(), 'g-')

		# fix the limits on plot
		xmin, xmax = fig.ax.get_xlim()
		ymin, ymax = fig.ax.get_ylim()
		plt.axis([xmin, xmax, ymin, ymax])

		# axes for parameter sliders
		self.Axis = {

			# key : axis     xpos , ypos + dy       , xsize, ysize
			line  : plt.axes([0.13, 0.05 + (k+1) * 0.03, 0.65, 0.02], axisbg = 'white')
			for k, line in enumerate( self.Params['L1'].keys() )
		}

		# add an axis to make adjustments to the continuum
		self.Axis['Continuum'] = plt.axes([0.13, 0.05, 0.65, 0.02], axisbg='white')

		# create the sliders for the parameters
		self.Slider = {

			# Slider `key` and widget
			param : widgets.Slider(

				self.Axis[param],    # which axis to put slider on
				param,               # name of parameter (and the slider)
				self.Minimum(param), # set the minimum of the slider
				self.Maximum(param), # set the maximum of the slider
				valinit = self.Params['L1'][param] # initial value
			)

			# create a slider for each parameter
			for param in self.Params['L1'].keys()
		}

		# create the slider for the continuum (10% adjustments)
		self.Slider['Continuum'] = widgets.Slider(self.Axis['Continuum'], 'Continuum',
			0.90 * self.continuum_level, 1.10 * self.continuum_level,
			valinit = self.continuum_level)

		# connect sliders to update function
		for slider in self.Slider.values():
			slider.on_changed(self.Update)

		# create axis for radio buttons
		self.RadioAxis = plt.axes([0.85, 0.01, 0.1, 0.2],
			axisbg = 'white', frameon = False)

		# create the radio button widget
		self.Radio = widgets.RadioButtons(self.RadioAxis,
			tuple(['L' + str(i+1) for i in range(self.numlines)]), active = 0)

		# connect the radio button to it's update function
		self.Radio.on_clicked(self.ToggleComponent)

		# set current component as the first line
		self.current_component = 'L1'

		# make currently selected line bold
		self.Component['L1'].set_linewidth(2)

		# variable allows for the radio buttons to not screw-up the sliders/graphs
		# when switching between lines
		self.stall = False

	def Parameterize(self, function, loc):
		"""
		Choose the initial parameters of `function` given the peak `loc`.
		"""

		if function == 'Voigt':

			return {

					# from the Lorentzian
					'Gamma': 0.50 * self.domainsize / self.numlines,

					# from the Gaussian
					'Sigma': 0.25 * self.domainsize / self.numlines,

					# center of the line
					'Peak' : loc[0],

					# depth of the line
					'Depth': self.continuum[ loc[0] ].value - loc[1]

				}

		if function == 'Lorentzian':

			return {

					# FWHM of the Lorentzian
					'FWHM': 0.5 * self.domainsize / self.numlines,

					# center of the line
					'Peak': loc[0],

					# depth of the line
					'Depth': self.continuum[ loc[0] ].value - loc[1]
				}

		elif function == 'Gaussian':

			return {

					# FWHM / Ln(2) of Gaussian
					'FWHM': 0.5 * self.domainsize / self.numlines,

					# center of line
					'Peak': loc[0],

					# depth of the line
					'Depth': self.continuum[ loc[0] ].value - loc[1]
				}

		else: raise ProfileError('From FittingGUI.Parameterize(), the only '
			'currently implemented functions are the `Lorentzian` and '
			'the `Gaussian`!')

	def SetFunction(self, function):
		"""
		Return which Evaluate result to use given profile function.
		"""

		options = {

			# local versions of functions handle parameters
			'Voigt'     : self.__Voigt,
			'Lorentzian': self.__Lorentzian,
			'Gaussian'  : self.__Gaussian
		}

		if function not in options:
			raise ProfileError('From FittingGUI.SetFunction(), the only '
			'currently implemented functions are the `Voigt`, `Lorentzian`, and '
			'the `Gaussian`!')

		return options[function]

	def __Gaussian(self, x, **params):
		"""
		A Gaussian profile. See ..Algorithms.Functions.Gaussian
		"""
		return Gaussian(

				x,

				# Amplitude of the Gaussian
				params['Depth'],

				# center
				params['Peak'],

				# sigma = FWHM / 2 sqrt{2 log 2}
				params['FWHM'] / 2.3548200450309493
			)


	def __Lorentzian(self, x, **params):
		"""
		A Lorentzian profile. See ..Algorithms.Functions.Lorentzian
		"""
		return params['Depth'] * Lorentzian(

				x,

				# location of the peak
				params['Peak'],

				# gamma *is* the FWHM
				params['FWHM']
			)

	def __Voigt(self, x, **params):
		"""
		The Voigt profile. See ..Algorithms.Functions.Voigt
		"""
		return Voigt(

				x,

				# depth of the line
				params['Depth'],

				# location of the peak
				params['Peak'],

				# sigma from the Gaussian
				params['Sigma'],

				# gamma from the Lorentzian
				params['Gamma']
			)

	def SuperPosition(self):
		"""
		Superposition of each line component
		"""
		# emtpy result
		combined = np.zeros(np.shape(self.x))

		# additional lines
		for parameters in self.Params.values():
			combined += self.Evaluate(self.x, **parameters)

		return combined

	def Minimum(self, param):
		"""
		Set the lower bound on the `param`eter for it's slider.
		"""
		if param == 'Peak':

			return self.x[0]

		elif param == 'FWHM':

			# don't divide by zero!!!!
			return 1e-6

		elif param == 'Depth':

			return 1e-6

		elif param == 'Sigma':

			return 1e-6

		elif param == 'Gamma':

			return 1e-6

		else: raise ProfileError('From FittingGUI.Minimum(), `{}` is not '
			'currently implemented as a parameter to set the minumum '
			'for!'.format(param))

	def Maximum(self, param):
		"""
		Set the upper bound on the `param`eter for it's slider.
		"""
		if param == 'Peak':

			return self.x[-1]

		elif param == 'FWHM':

			return 0.9 * self.domainsize

		elif param == 'Depth':

			return 1.5 * self.continuum.max()

		elif param == 'Sigma':

			return 0.9 * self.domainsize / self.numlines

		elif param == 'Gamma':

			return 0.9 * self.domainsize / self.numlines

		else: raise ProfileError('From FittingGUI.Maximum(), `{}` is not '
			'currently implemented as a parameter to set the maximum '
			'for!'.format(param))

	def Update(self, val):
		"""
		Cycle thru Sliders and update Parameter dictionary. Re-draw graphs.
		"""

		if not self.stall:

			# the currently selected line component
			line = self.current_component

			# update the appropriate parameters in the dictionary
			for parameter, slider in self.Slider.items():

				if parameter == 'Continuum':
					self.y = self.continuum + (slider.val - self.continuum_level)

				else:
					self.Params[line][parameter] = slider.val

			# update the appropriate graph data, based on new parameters
			for line, parameters in self.Params.items():
				self.Component[line].set_ydata(self.y - self.Evaluate(self.x, **parameters))

			# update the super-imposed graphs
			self.Combination.set_ydata(self.y - self.SuperPosition())

			# push updates to graph
			self.fig.canvas.draw_idle()

	def ToggleComponent(self, label):
		"""
		Toggle function for the radio buttons. Switch between line components
		`L1`, `L2`, etc. Update the sliders to reflect changing parameters.
		"""

		# reassign the current component that is selected
		self.current_component = label

		# don't update the parameter dictionary or draw the graph!!!!!!
		self.stall = True

		# make current feature bold and the rest thin
		for line in self.Component.keys():
			if line == label:
				self.Component[line].set_linewidth(2)
			else:
				self.Component[line].set_linewidth(1)

		# update the sliders to reflect the current component
		for parameter, slider in self.Slider.items():
			if parameter != 'Continuum':
				slider.set_val(self.Params[label][parameter])

		# give control back to sliders
		self.stall = False

		# push updates to graph
		self.fig.canvas.draw_idle()

	def GetSpectra(self):
		"""
		Return the current graphs as Spectrum objects.
		"""
		return [

				Spectrum(
						# build a spectrum based on numpy arrays from plot
						gui.Component[graph].get_ydata() * gui.line.data.unit,
						gui.x * gui.line.wave.unit
					)

				# for all `L1`, `L2`, etc...
				for graph in sorted( gui.Component.keys() )
			]

	def GetContinuum(self):
		"""
		Return the current graph of the continuum.
		"""
		return Spectrum(

				# data from the graph, with applied units
				self.y * self.line.data.unit,
				self.x * self.line.wave.unit
			)

def MultiFit(splot, measurements=True, resolution=1e5, function='Voigt', **kwargs):
	"""
	The MultiFit routine takes a `splot` figure (type SPlot) and allows the
	user to interactively fit line profiles. `splot` may optionally be of type
	Spectrum, in which case a SPlot figure will be created for you. This function
	creates a *FittingGUI* object (not documented here) which uses the *Profile.Extract()*
	routine first (*kwargs* as passed to this function). As in the Extract routine, the user
	will select four points that mark the boundaries of the line (blended or otherwise)
	and some surrounding continuum to sample. The *KernelFit1D* routine is applied with the
	user specified *bandwidth* to smooth the noise out of the continuum and interpolate, of
	type *kind*, across the line gap. After this, the user is asked to further select points
	marking the peaks of the line(s). The number of selection points indicates the number
	of lines to be fit. If you wish to deblend two or more lines, select all suspected
	locations. These are not only used to determine the number of lines to fit but to make
	initial guesses at the parameters for each line and determine reasonable scales for the
	sliders.

	After these selections, a slider for each parameter appears along with radio buttons for
	each line. The user can interactively adjust the parameter(s) for a given line by
	selecting it (the radio buttons) and dragging a slider. Each line is represented by a
	black, dashed line in the plot. The currently selected line is bold. The combined
	(i.e., blended) lines are plotted with a solid green line.

	Each slider is labeled on the left side with the name of the parameter it controls. At
	the right of each slider is the current value of that parameter. The radio buttons are
	labeled "L1", "L2", etc. for each line.

	This routine returns a list of Spectrum objects. The length of the list is equal to the
	number of lines that were fit plus one for the continuum. This routine functions as both
	a single line fitting tool as well as a deblending tool. By default, the final
	parameterizations are attached to each spectrum as a dictionary. Also by default, various
	measurements from the fitted profiles are attached to each returned Spectrum object. These
	include the equivelent width (attached as `.ew`) and column depth (attached as `N`). This
	behavior can be suppressed by giving the keyword argument `measurements = False`. The
	measurements are of type `Measurement` (..Framework.Measurement.Measurement). This is
	so the uncertainty for each measurement can be attached to its result.

	Example:
	```
	Lines = Profile.MultiFit(fig, bandwidth = u.Angstrom / 10)
	Lines[1].N
	`Equivelent Width`: (1.4567 ^+ 0.0123 _- 0.0213) x 10^-3 Angstrom

	```
	"""

	# running the user interface
	gui = FittingGUI(splot, function=function, **kwargs)

	input('\n Press <Return> when you are finished fitting lines ...')

	# attach the parameter dictionaries to each spectrum
	lines = gui.GetSpectra()
	for a, parameterization in enumerate( sorted(gui.Params.keys()) ):
		lines[a].parameters = gui.Params[parameterization]

	if measurements:

		print('\n Computing the Equivalent Widths ... ', end = '')

		continuum = lines[-1]

		# the ratio of the line to the continuum, percent absorption (inverted I guess)
		I0_I = [ continuum / spectrum for spectrum in lines[:-1] ]

		# RMS error from the continuum fit (i.e., KernelFit1D in Profile.Extract())
		rms = gui.rms

		# the error in each line is from both the continuum noise and the counts
		line_errors = [

				# the error is added in quadrature,
				# error in the continuum is rms
				# error in the line is rms * sqrt(I / I0)
				Spectrum( rms * np.sqrt(1 + 1 / absorption.data) * continuum.data.unit,
					continuum.wave)
			]

		# the Equivelent Width
		EWs = [
				simps(spectrum.wave, continuum - 1 / absorption) * spectrum.wave.unit
				for spectrum, absorption in zip(lines, I0_I)
			]
