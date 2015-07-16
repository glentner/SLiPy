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
from scipy.integrate import simps as Integrate

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import widgets

from astropy import units as u
from astropy.constants import m_e, e, c

from .. import SlipyError

from .Observatory import Observatory
from .Spectrum import Spectrum, SpectrumError
from .Plot import SPlot, PlotError

from ..Framework.Options import Options, OptionsError
from ..Framework.Measurement import Measurement
from ..Algorithms.Functions import Gaussian, Lorentzian, Voigt, InvertedLorentzian
from ..Algorithms.KernelFit import KernelFit1D
from ..Data import Atomic

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
            marker = 'o', s=75, edgecolor='red', facecolor='None', lw=2)
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
    Given `splot` of type SPlot, the user selects four points on the
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
        'the Profile.AutoFit() routine!')

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
            'Profile.AutoFit()!')

        try:
            for a, parameter in enumerate(params):
                if type(parameter) is type(InvertedLorentzian):
                    # replace parameter with function evaluation
                    params[a] = parameter(x_inner, y_inner)

        except TypeError as err:
            print(' --> TypeError:', err)
            raise ProfileError('Profile.AutoFit() failed to call user functions '
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
			'bandwidth' : 0.1*u.nm  # user should provide this!
			# 'rms'       : False     # measure the RMS of the line, continuum
			})

		kind      = options('kind')
		bandwidth = options('bandwidth')
		# rms       = options('rms')

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

	# if not rms:
	#     return line, continuum

	continuum.rms  = np.sqrt( np.sum( (cont_outside * yc.unit - yc)**2 ) / len(yc) )
	continuum.rms *= continuum.data.unit

	return line, continuum


def OpticalDepth(line, continuum, error=None, boost=None):
	"""
	Given an absorption `line` and its background `continuum`, compute the
	apparent optical depth Spectrum. Both the line and the continuum must
	be of Spectrum type. The continuum will be resampled to the pixel space of the
	line if it is not currently.

	If provided an `error` spectrum, it should either have the same units as the `line`
	or be in percent units (they should have identical wavelength arrays). An upper and
	lower uncertainty will be computed by adding and subtracting before the calculation.
	Further, if an `rms` value is attached to the `continuum` giving the percent error in
	the continuum fit, this will be added in quadrature to the error spectrum before hand.

	`boost` allows for artificially increasing the resolution by resampling to more pixels.
	If the spectrum is not of sufficiently high resolution, the integration could suffer
	numerical errors. If provided, this should be a percentage to increase the resolution
	(e.g., `boost=2` would double the resolution by interpolating in between the pixels).

	This function returns a Spectrum with `upper` and `lower` bounds attached.
	"""
	if not isinstance(line, Spectrum):
		raise ProfileError('OpticalDepth() expects the first argument to be of '
		'`Spectrum` type!')

	if not isinstance(continuum, Spectrum):
		raise ProfileError('OpticalDepth() expects the second argument to be of '
		'`Spectrum` type!')

	line      = line.copy()
	continuum = continuum.copy()

	if error:
		if not isinstance(error, Spectrum):
			raise ProfileError('OpticalDepth() expects the `error` to be of `Spectrum` type!')
		if error.data.unit not in [line.data.unit, u.percent]:
			raise ProfileError('OpticalDepth() expects the `error` spectrum to have either '
			'the same units as the `line` or units of `percent`.')

		error = error.copy()

		if error.data.unit == u.percent:
			if np.logical_and(error.data.value < 0, error.data.value > 100).any():
				raise ProfileError('OpticalDepth() was given an `error` spectrum in '
				'units of `percent` with values outside of 0 and 100!')

			error.resample(line)
			error.data = (error.data * line.data).to(line.data.unit)

	if hasattr(continuum, 'rms'):
		if not hasattr(continuum.rms, 'unit') or (continuum.rms.unit not in
			[continuum.data.unit, u.percent]):
			raise ProfileError('OpticalDepth() expects a `continuum` with an `rms` '
			'attribute to have the same units or `percent` units.')

		rms = continuum.rms
		if rms.unit == u.percent:
			rms = rms * continuum.data
			rms = rms.to(continuum.data.unit)

		if not error:
			error = rms

		else:
			# add the two error components in quadrature
			error.resample(line)
			error = Spectrum(np.sqrt(rms**2 + error.data**2) * line.data.unit, line.wave)

	# boost the resolution of the specrum if requested
	if boost: line.resample( line.wave[0], line.wave[-1], len(line) * boost )

	# resample the continuum spectrum, nothing happens if they are already the same
	continuum.resample(line)

	# compute the apparent optical depth
	tau = Spectrum(np.log((continuum / line).data.decompose()), line.wave)

	if error:
		tau.lower = Spectrum(np.log((continuum / (line - error)).data.decompose()), line.wave)
		tau.upper = Spectrum(np.log((continuum / (line + error)).data.decompose()), line.wave)

	return tau

def EquivalentWidth(line, continuum, error=None, boost=None):
	"""
	Given an absorption `line` and its background `continuum`, compute the
	equivalent width, `W`, of the feature. Both the line and the continuum must
	be of Spectrum type. The continuum will be resampled to the pixel space of the
	line if it is not currently.

	If provided an `error` spectrum, it should either have the same units as the `line`
	or be in percent units (they should have identical wavelength arrays). An upper and
	lower uncertainty will be computed by adding and subtracting before the calculation.
	Further, if an `rms` value is attached to the `continuum` giving the percent error in
	the continuum fit, this will be added in quadrature to the error spectrum before hand.

	`boost` allows for artificially increasing the resolution by resampling to more pixels.
	If the spectrum is not of sufficiently high resolution, the integration could suffer
	numerical errors. If provided, this should be a percentage to increase the resolution
	(e.g., `boost=2` would double the resolution by interpolating in between the pixels).

	Integration is performed using the composite Simpson`s rule (scipy.integrate.simps)
	This function returns a `Measurement` object (..Framework.Measurement.Measurement).
	"""
	tau = OpticalDepth(line, continuum, error=error, boost=boost)
	W   = Integrate(1 - np.exp(-tau.data.decompose()), tau.wave) * line.wave.unit

	if error:

		upper = Integrate(1 - np.exp(-tau.lower.data.decompose()), tau.wave) * line.wave.unit
		lower = Integrate(1 - np.exp(-tau.upper.data.decompose()), tau.wave) * line.wave.unit

		# join upper/lower into +/- set
		uncertainty = np.array([(upper - W).value, (lower - W).value]) * tau.wave.unit

	else: uncertainty = None

	return Measurement(W, error = uncertainty, name = 'Equivalent Width',
		notes = 'Measured using Profile.EquivalentWidth() from SLiPy')


def ColumnDensity(line, continuum, ion, error=None, boost=None, weakline=None, integrate=True):
	"""
	Given an absorption `line` and its background `continuum`, compute the
	apparent column density, `N`, of the feature. Both the line and the continuum must
	be of Spectrum type. The continuum will be resampled to the pixel space of the
	line if it is not currently. An `ion` must be provided (type Atomic.Ion) with
	members, `fvalue` (oscillator strength) and `wavelength` (will be converted to
	Angstroms).

	If provided an `error` spectrum, it should either have the same units as the `line`
	or be in percent units (they should have identical wavelength arrays). An upper and
	lower uncertainty will be computed by adding and subtracting before the calculation.
	Further, if an `rms` value is attached to the `continuum` giving the percent error in
	the continuum fit, this will be added in quadrature to the error spectrum before hand.

	With `integrate` set to True (the default) the integrated apparent column density will
	be returned (as type Measurement and in units of cm-2). Otherwise, a Spectrum is
	returned as a function of wavelength and the `.upper` and `.lower` bounds are attached.

	`boost` allows for artificially increasing the resolution by resampling to more pixels.
	If the spectrum is not of sufficiently high resolution, the integration could suffer
	numerical errors. If provided, this should be a percentage to increase the resolution
	(e.g., `boost=2` would double the resolution by interpolating in between the pixels).

	With a `weakline` provided (the results of this function on the weaker absorption
	line for a doublet) a correction gives an approximate `true` column density (attached
	to the returned Measurement as `.true`), see Savage & Sembach, 1991.

	Integration is performed using the composite Simpson`s rule (scipy.integrate.simps)
	This function returns a `Measurement` object (..Framework.Measurement.Measurement).
	"""
	if (not isinstance(ion, Atomic.Ion) or not hasattr(ion, 'fvalue') or
		not hasattr(ion, 'wavelength')): raise ProfileError("From ColumnDensity(), `ion` is "
		"expected to be of type Atomic.Ion and have members `fvalue` and `wavelength`!")

	if weakline and not hasattr(weakline, 'unit'):
		raise ProfileError("From ColumnDensity(), `weakline` is expected to have units!")

	const   = 1 / (ion.fvalue * ion.wavelength.to(u.AA)**2 * np.pi * e.emu**2 / (m_e.si *
				c.to(u.km / u.second)**2))
	tau     = OpticalDepth(line, continuum, error=error, boost=boost)
	N       = const.value * tau       / (u.cm**2 / tau.wave.unit)
	N.upper = const.value * tau.upper / (u.cm**2 / tau.wave.unit)
	N.lower = const.value * tau.upper / (u.cm**2 / tau.wave.unit)

	if not integrate and not weakline:
		return N

	elif weakline and not integrate:
		raise ProfileError("From ColumnDensity(), you gave `integrate` as False but we "
		"need to integrate to solve for the correct N using the `weakline`!")

	upper = Integrate(N.upper.data, N.wave) / u.cm**2
	lower = Integrate(N.lower.data, N.wave) / u.cm**2
	N     = Integrate(N.data, N.wave) / u.cm**2
	uncertainty = np.array([(N - upper).value, (lower - N).value]) / u.cm**2

	N = Measurement(N, name="Column Density (Apparent)",
		error=uncertainty, notes="Measured using Profile.ColumnDensity() from SLiPy")

	if not weakline:
		return N

	# attach the `correction` given the already computed `weakline` column density.
	diff = np.log10(weakline.to(1 / u.cm**2).value) - np.log10(N.value)
	if diff < 0.0 or diff > 0.24:
		raise ProfileError("From ColumnDensity(), the difference between the doublet "
		"lines is too great to solve for a correction using published results!")

	# Table 4: Column Density Corrections for an Isolated Gaussian Component
	# Savage & Sembach (1991)
	table = np.array([  [0.000, 0.000], [0.010, 0.010], [0.020, 0.020], [0.030, 0.030],
		[0.040, 0.040], [0.050, 0.051], [0.060, 0.061], [0.070, 0.073], [0.080, 0.085],
		[0.090, 0.097], [0.100, 0.111], [0.110, 0.125], [0.120, 0.140], [0.130, 0.157],
		[0.140, 0.175], [0.150, 0.195], [0.160, 0.217], [0.170, 0.243], [0.180, 0.273],
		[0.190, 0.307], [0.200, 0.348], [0.210, 0.396], [0.220, 0.453], [0.230, 0.520],
		[0.240, 0.600]])

	N.correction = interp1d(table[:,0], table[:,1], kind='cubic')(diff)
	N.true       = 10**(np.log10(weakline.to(1/u.cm**2).value) + N.correction) / u.cm**2
	return N


class FittingGUI:
	"""
	Graphical Interface (keeps references alive) for fitting analytical
	profiles to spectral data.
	"""
	def __init__(self, splot, obs=None, ion=None, function='Voigt', **kwargs):
		"""
		Build widget elements based on a spectrum plotted in `splot` (type SPlot).
		`splot` may also be a Spectrum object (for which a SPlot will be created).

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
			raise ProfileError('Unrecognized option given to FittingGUI.__init__()!')

		if function not in ['Voigt', 'Lorentzian', 'Gaussian']:
			raise ProfileError('The only currently implemented functions '
			'for fitting profiles are `Voigt`, `Lorentzian` and `Gaussian`!')

		# grab and/or create the SPlot
		if isinstance(splot, Spectrum):
			splot = SPlot(splot, marker='k-', label='spectrum')

		elif not isinstance(splot, SPlot):
			raise ProfileError('FittingGUI() expects the `splot` argument to '
			'either be a Spectrum or a SPlot!')

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# if the profile `function` is `Voigt` and we have necessary parameters,
		# show a preview for `b` and `N` and `v`.
		self.any_line_parameters = True if obs or  ion else False
		self.has_line_parameters = True if obs and ion else False
		if self.any_line_parameters and not self.has_line_parameters:
			raise ProfileError('From FittingGUI.__init__(), if the absoption line '
			'parameters `obs` or `ion` are provided, both must be given!')

		if self.has_line_parameters:

			if function != 'Voigt':
				raise ProfileError('From FittingGUI.__init__(), in order to compute the '
				'broadening from the instrument profile of the `obs`ervatory and the '
				'column density for the `ion`, it is expected that we try to fit a `Voigt` '
				'profile `function`!')

			if not isinstance(obs, Observatory):
				raise ProfileError('From FittingGUI.__init__(), if `obs` is provided it should '
				'be of type Observatory!')

			if not hasattr(obs, 'resolution'):
				raise ProfileError('From FittingGUI.__init__(), if an observatory `obs` is '
				'provided, it needs to have a member `resolution`!')

			if not hasattr(ion, 'wavelength') or not hasattr(ion.wavelength, 'unit'):
				raise ProfileError('From FittingGUI.__init__(), the provided `ion` either '
				'did not have a `wavelength` attribute or it has one without units!')

			if not hasattr(ion, 'fvalue'):
				raise ProfileError('From FittingGUI.__init__(), the provide `ion` does not '
				'have an `fvalue` (oscillator strength) attribute!')

			if hasattr(ion.fvalue, 'unit') and ion.fvalue.unit != u.dimensionless_unscaled:
				raise ProfileError('From FittingGUI.__init__(), the osciallator strength, '
				'`fvalue` for an ion must be a dimensionless Quantity!')

			if not hasattr(ion, 'A') or not ion.A:
				raise ProfileError('From FittingGUI.__init__(), the provided `ion` does not '
				'have an `A` (transition probability) attribute!')

			if hasattr(ion.A, 'unit') and ion.A.unit != u.Unit('s-1'):
				raise ProfileError('From FittingGUI.__init__(), the transition probability, '
				'`A`, from the `ion` must be in units of `s-1` if units are present!')

			# the instrument broadening is from R = lambda / delta_lambda
			# FWHM = 2 sqrt( 2 log 2) sigma for the Gaussian instrument profile
			self.R = obs.resolution
			self.sigma_instrument = (ion.wavelength / self.R) / (2 * np.sqrt(2 * np.log(2)))
			self.sigma_instrument_squared = self.sigma_instrument.value ** 2

			# save parameters for later
			self.wavelength = ion.wavelength
			self.fvalue     = ion.fvalue
			self.A          = ion.A

			# the FWHM of the intrinsic line profile (Lorentzian) is proportional to the
			# transition probability (Einstein coeff.) `A`...
			# convert from km s-1 to wavelength units
			self.gamma = (ion.wavelength * (ion.wavelength * ion.A / (2 * np.pi)).to(u.km / u.s) /
				c.si).to(ion.wavelength.unit).value

			# the leading constant in the computation of `N` (per Angstrom per cm-2)
			# self.leading_constant = (m_e.si * c.si / (np.sqrt(np.pi) * e.si**2 * ion.fvalue *
			# 	ion.wavelength.to(u.Angstrom))).value
			self.leading_constant = 1 / (ion.fvalue * ion.wavelength.to(u.AA)**2 * np.pi *
				e.emu**2 / (m_e.si * c.to(u.km/u.s)**2)).value

			# setting `function` to `ModifiedVoigt` only makes a change in how the `Voigt`
			# profile is evaluated by using `self.gamma` instead of self.Params[...]['Gamma']
			function = 'ModifiedVoigt'

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		print('\n We need to extract the lines from the continuum before we '
		'begin the fitting process.')

		# extract the line, continuum, rms from the spectrum
		self.line, self.continuum = Extract(splot, bandwidth=bandwidth, kind=kind)
		self.rms = self.continuum.rms

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

		# final spectral line profile is convolution of the LSP and the
		# line profile function. Gaussian on Gaussian, Gaussian on Lorentzian,
		# etc ... Set the functional form given requested profile function
		self.Evaluate = self.SetFunction(function)

		# initial guesses for parameters given the line locations,
		# containing `L1`, `L2`, etc ... the values of which are themselves
		# dictionaries of parameters (e.g., `FWHM`, `Depth`, etc...) whose values
		# are the initial guesses for those parameters given the location of
		# it`s peak and the number of peaks within the `line`s domain.
		self.Params = {

				'L' + str(line + 1) : self.Parameterize(function, loc)
				for line, loc in enumerate(points[:, 2:-2].T)

			}

		# grab the actual Figure object and it's axis to keep references
		self.splot = splot
		self.fig   = splot.fig
		self.ax    = splot.ax

		# refresh image, but keep any changes in axis limits
		splot.xlim( *splot.ax.get_xlim() )
		splot.ylim( *splot.ax.get_ylim() )
		splot.draw()

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

		# add plot of the continuum
		self.Continuum, = plt.plot(self.x, self.y, 'r-', lw=2)

		# add plot of superposition of each component
		self.Combination, = plt.plot(self.x, self.y - self.SuperPosition(), 'g-')

		# fix the limits on plot
		xmin, xmax = splot.ax.get_xlim()
		ymin, ymax = splot.ax.get_ylim()
		plt.axis([xmin, xmax, ymin, ymax])

		# axes for parameter sliders
		self.Axis = {

			# key : axis     xpos , ypos + dy       , xsize, ysize
			line  : plt.axes([0.10, 0.05 + (k+1) * 0.03, 0.65, 0.02], axisbg = 'white')
			for k, line in enumerate( self.Params['L1'].keys() )
		}

		# add an axis to make adjustments to the continuum
		self.Axis['Continuum'] = plt.axes([0.10, 0.05, 0.65, 0.02], axisbg='white')

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

		# add the initial text for `N` and `b` if applicable.
		# text is along 20% below the y-axis
		if self.has_line_parameters:
			# self.b, self.N, self.v = self.solve_b(), self.solve_N(), self.solve_v()
			self.preview = self.ax.text(xmin, ymin - 0.1 * (ymax - ymin),
				'v: {:>10.4f}        km s-1            b: {:>10.4f}    km s-1'
				'       W: {:>10.4f}\nN:   {:>10.4e}    cm-2      tau_0: {:>10.4f}'.format(
				self.solve_v().value, self.solve_b().value, self.solve_W(), self.solve_N().value,
				self.solve_tau0()), va = 'top')

		# display the text
		self.fig.canvas.draw_idle()

	def Parameterize(self, function, loc):
		"""
		Choose the initial parameters of `function` given the peak `loc`.
		"""
		if function == 'Voigt':
			return {'Gamma': 0.10 * self.domainsize / self.numlines,
				'Sigma': 0.20 * self.domainsize / self.numlines, 'Peak' : loc[0],
				'Depth': self.continuum[ loc[0] ].value - loc[1]}

		if function == 'Lorentzian':
			return {'FWHM': 0.5 * self.domainsize / self.numlines, 'Peak': loc[0],
				'Depth': self.continuum[ loc[0] ].value - loc[1]}

		elif function == 'Gaussian':
			return {'FWHM': 0.5 * self.domainsize / self.numlines, 'Peak': loc[0],
				'Depth': self.continuum[ loc[0] ].value - loc[1]}

		elif function == 'ModifiedVoigt':
			return {'Sigma': (self.Maximum('Sigma') - self.Minimum('Sigma'))/2,
				'Peak' : loc[0], 'Depth': self.continuum[ loc[0] ].value - loc[1]}

		else:
			raise ProfileError('From FittingGUI.Parameterize(), the only '
			'currently implemented functions are the `Gaussian`, `Lorentzian`, and '
			'`Voigt` profiles!')

	def SetFunction(self, function):
		"""
		Return how to `Evaluate` the profile given the `function`.
		"""
		options = {'Voigt': self.__Voigt, 'Lorentzian': self.__Lorentzian,
			'Gaussian': self.__Gaussian, 'ModifiedVoigt': self.__ModifiedVoigt}

		if function not in options:
			raise ProfileError('From FittingGUI.SetFunction(), the only '
			'currently implemented functions are the `Voigt`, `Lorentzian`, and '
			'the `Gaussian` profiles!!')

		return options[function]

	def __Gaussian(self, x, **params):
		"""
		A Gaussian profile. See ..Algorithms.Functions.Gaussian
		"""
		return Gaussian(x, params['Depth'], params['Peak'], params['FWHM'] / 2.3548200450309493)

	def __Lorentzian(self, x, **params):
		"""
		A Lorentzian profile. See ..Algorithms.Functions.Lorentzian
		"""
		return params['Depth'] * Lorentzian(x, params['Peak'], params['FWHM'])

	def __Voigt(self, x, **params):
		"""
		The Voigt profile. See ..Algorithms.Functions.Voigt
		"""
		return Voigt(x, params['Depth'], params['Peak'], params['Sigma'], params['Gamma'])

	def __ModifiedVoigt(self, x, **params):
		"""
		This is the Voigt profile, but `gamma` was already set.
		"""
		return Voigt(x, params['Depth'], params['Peak'], params['Sigma'], self.gamma)

	def SuperPosition(self):
		"""
		Superposition of each line component (blended line)
		"""
		combined = np.zeros(np.shape(self.x))
		for parameters in self.Params.values():
			combined += self.Evaluate(self.x, **parameters)
		return combined

	def Minimum(self, param):
		"""
		Set the lower bound on the `param`eter for its slider.
		"""
		if param == 'Peak':
			return self.x[0]

		elif param == 'FWHM':
			return 1e-6

		elif param == 'Depth':
			return 1e-6

		elif param == 'Sigma':
			if self.Evaluate != self.__ModifiedVoigt:
				return 1e-6
			else:
				return self.sigma_instrument.value

		elif param == 'Gamma':
			return 1e-6

		else:
			raise ProfileError('From FittingGUI.Minimum(), `{}` is not '
			'currently implemented as a parameter to set the minumum '
			'for!'.format(param))

	def Maximum(self, param):
		"""
		Set the upper bound on the `param`eter for its slider.
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

		else:
			raise ProfileError('From FittingGUI.Maximum(), `{}` is not '
			'currently implemented as a parameter to set the maximum '
			'for!'.format(param))

	def Update(self, val):
		"""
		Cycle thru Sliders and update Parameter dictionary. Re-draw graphs.
		"""
		# `self.stall` suppresses this procedure while inside `ToggleComponent`
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

			# update the continuum graph
			self.Continuum.set_ydata(self.y)

			# update the combined graph
			self.Combination.set_ydata(self.y - self.SuperPosition())

			# update the displayed `N` and `b` if present
			if self.has_line_parameters:
				self.Update_Preview()

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

		# update the displayed `b`, `N` and `z` if present
		if self.has_line_parameters:
			self.Update_Preview()

		# give control back to sliders
		self.stall = False

		# push updates to graph
		self.fig.canvas.draw_idle()

	def solve_b(self):
		"""
		Given the current line component, solve for the broadening parameter, `b`,
		given the instrument profile and the observed broadening.
		"""
		# b = sqrt(2) * sigma_v
		sigma_obs = self.Params[self.current_component]['Sigma']
		sigma_v   = np.sqrt(sigma_obs**2 - self.sigma_instrument_squared) * self.wavelength.unit
		b = np.sqrt(2) * (c * sigma_v / self.wavelength)
		return b.to(u.km / u.second)
		# return (c.si * (1.4142135623730951 *
		# 	np.sqrt(self.Params[self.current_component]['Sigma']**2
		# 	- self.sigma_instrument_squared) * self.wavelength.unit) / self.wavelength).to(u.km /
		# 	u.second)

	def solve_tau0(self):
		"""
		Solve for the apparent optical depth at line center.
		"""
		line_data = self.Component[self.current_component].get_ydata()
		line_wave = self.Component[self.current_component].get_xdata()
		line_center = line_data.argmin()

		return np.log(self.y[line_center] / line_data[line_center])

	def solve_W(self):
		"""
		Solve for the equivalent width of the currently selected line.
		"""
		line_data = self.Component[self.current_component].get_ydata()
		line_wave = self.Component[self.current_component].get_xdata()

		# apparent optical depth
		tau = np.log(self.y / line_data)
		return Integrate(1 - np.exp(-tau), line_wave) * self.wavelength.unit

	def solve_N(self):
		"""
		Solve for the column density of the currently selected line component.
		"""
		# equivalent width (dimensionless)
		W = self.solve_W() / self.wavelength

		# m_e c^2 / pi e^2 \approx 1.13e12 cm-1 ???? How is that?
		const = (1.13e12 / u.cm) / (self.fvalue * self.wavelength.to(u.cm))
		return const * W

	def solve_v(self):
		"""
		Solve for the velocity of the line given the expected wavelength.
		The result is returned in km s-1.
		"""
		line_center = self.Params[self.current_component]['Peak'] * self.wavelength.unit
		return c.to(u.km / u.second) * (line_center - self.wavelength) / self.wavelength

	def Update_Preview(self):
		"""
		Re-compute the `b`, `N`, and `v` values, update the text in the plot.
		"""
		# self.b, self.N, self.v = self.solve_b(), self.solve_N(), self.solve_v()
		self.preview.set_text('v: {:>10.4f}        km s-1            b: {:>10.4f}    km s-1'
		'       W: {:>10.4f}\nN:   {:>10.4e}    cm-2      tau_0: {:>10.4f}'.format(
		self.solve_v().value, self.solve_b().value, self.solve_W(), self.solve_N().value,
		self.solve_tau0()))

	def GetSpectra(self):
		"""
		Return the current graphs as Spectrum objects.
		"""
		return [Spectrum(self.Component[line].get_ydata() * self.line.data.unit,
			self.x * self.line.wave.unit) for line in sorted(self.Component.keys())]

	def GetContinuum(self):
		"""
		Return the current graph of the continuum.
		"""
		continuum = Spectrum(self.y * self.line.data.unit, self.x * self.line.wave.unit)
		continuum.rms = self.rms
		return continuum

	def kill(self):
		"""
		Close the plot to destroy widgets and restore the `splot` to its original state.
		"""
		plt.close(self.fig)
		self.splot.fig = plt.figure("Spectral-Plot (SLiPy)")
		self.splot.ax  = self.splot.fig.add_subplot(111)
		self.splot.draw()
		#del(self)

def MultiFit(splot, error=None, obs=None, ion=None, function='Voigt', measure=True,
	boost=None, **kwargs):
	"""
	The MultiFit routine takes a `splot` figure (type SPlot) and allows the
	user to interactively fit line profiles. `splot` may optionally be of type
	Spectrum, in which case a SPlot figure will be created for you. This function
	creates a *FittingGUI* object (not documented here) which uses the Profile.Extract()
	routine first (`kwargs` are passed to this function). As in the Extract routine, the user
	will select four points that mark the boundaries of the line (blended or otherwise)
	and some surrounding continuum to sample. The KernelFit1D (..Algorithms.KernelFit) routine
	is applied with the user specified `bandwidth` to smooth the noise out of the continuum
	and interpolate, of type `kind`, across the line gap. After this, the user is asked to
	further select points marking the peaks of the line(s). The number of selection points
	indicates the number of lines to be fit. If you wish to deblend two or more lines, select
	all suspected locations. These are not only used to determine the number of lines to fit
	but to make initial guesses at the parameters for each line and determine reasonable scales
	for the sliders.

	After these selections, a slider for each parameter appears along with radio buttons for
	each line. The user can interactively adjusts the parameter(s) for a given line by
	selecting it (the radio buttons) and dragging a slider. The parameters available to
	adjust depend on the `function` argument. There are currently three available line shapes:
	'Gaussian', 'Lorentzian', and 'Voigt'. See ..Algorithms.Functions for information on
	these.

	Each line is represented by a black, dashed curve in the plot. The currently selected line
	is bold. The combined (i.e., blended) line is plotted as a solid green curve.

	If an Observatory with a `.resolution` is provided via the `obs` argument then the
	thermal broadening parameter `b` can be computed (and displayed). This is only applicable
	with either a `Gaussian` or `Voigt` profile. Further, an `ion` needs to be provided in
	this scenario (type ..Data.Atomic.Ion) with an oscillator strength (fvalue), transition
	probability (A) and wavelength as members. With these data, we can interactively show
	the broadening prameter, the velocity, and the apparent column density during the fitting
	process. Note: the `gamma` slider will disappear in this context because it has already
	been determined via the transition probability.

	Each slider is labeled on the left side with the name of the parameter it controls. At
	the right of each slider is the current value of that parameter. The radio buttons are
	labeled "L1", "L2", etc. for each line.

	This routine returns both a list of Spectrum `lines` as well as the `continuum` that was
	fit using the FittingGUI. This functions acts as both a single line fitting tool as well
	as a deblending tool. By default, the final parameterizations are attached to each spectrum
	as a dictionary, `.parameters`.

	If `measure` is passed as True (the default), the equivalent width is computed and attached
	to each spectrum (`L1`, `L2`, etc...) and can be accessed with `.W`. If `b` and/or `N` are
	available for each line these will be attached as well, accessible as `.b` and `.N`
	respectively. This functionality can be suppressed by setting `measure` to False. `boost`
	and `error` are passed to EquivalentWidth() and ColumnDensity().
	"""

	# running the user interface
	gui = FittingGUI(splot, obs=obs, ion=ion, function=function, **kwargs)
	input('\n Press <Return> when you are finished fitting lines ...')

	# attach the parameter dictionaries to each spectrum
	lines = gui.GetSpectra()
	for a, parameterization in enumerate( sorted(gui.Params.keys()) ):
		lines[a].parameters = gui.Params[parameterization]

	# attach the continuum to the line list
	continuum = gui.GetContinuum()

	if not measure:
		return lines, continuum

	if gui.has_line_parameters:
		for line, key in zip(lines, sorted(gui.Params.keys())):
			# set the line to `L1`, `L2`, etc...
			gui.current_componenet = key
			# attach the measured quantities
			notes  = 'Measurement made using Profile.MultiFit() from SLiPy'

			line.W = EquivalentWidth(line, continuum, error=error, boost=boost)
			line.N = ColumnDensity(line, continuum, ion, error=error, boost=boost)
			line.b = Measurement(gui.solve_b(), name='Doppler Broadening Paramater', notes=notes)
			line.v = Measurement(gui.solve_v(), name='Velocity', notes=notes)

	gui.kill()

	return lines, continuum
