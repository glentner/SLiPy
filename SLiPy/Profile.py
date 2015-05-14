# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/SLiPy/Profile.py 

"""
Profile fitting tasks for spectra.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate.interpolate import interp1d

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from astropy import units as u

from .. import SlipyError
from .DataType import Spectrum, DataTypeError
from .Plot import SPlot, PlotError
from ..Framework.Options import Options, OptionsError

from ..Algorithms.Functions import Gaussian
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
    Select points from the `splot`. If splot should be of type SPlot 
    (or it can optionally be a Spectrum type, for which a SPlot will be
    created). The splot will be rendered and the user clicks on the 
    figure. When finished, return to the terminal prompt.
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

def Fit(splot):
	"""
	Given `splot` of type Plot.SPlot, the user selects two points on the 
    spectrum and
	"""
	try:

		# function parameter defaults
		options = Options( kwargs,
			{
				'pad': 0 # expantion of domain for evaluation of the profile function.
			})

		# function parameter assignments
		pad = options('pad')

		# check input arguments
		if type(spectrum) is not Spectrum:
			raise ProfileError('Fit expects type Spectrum.')
		if xmin < spectrum.wave[0] or xmax > spectrum.wave[-1]:
			raise ProfileError('Out of domain for `spectrum`.')

		# make local copy of spectrum
		local = spectrum.copy()

		# extract relavent data from spectrum
		wave, data = local.wave, local.data
		data = data[ wave[ wave < xmax ] > xmin ]
		wave = wave[ wave[ wave < xmax ] > xmin ]

		# initial guess of parameters
		p0 = [ 1 - data.min(), wave[ data.argmin() ], data.std() ]

		# fit curve
		coeff, var_matrix = curve_fit(igauss, wave, data, p0=p0)

		# create Spectrum of profile
		x = local.wave
		x = x[ x[ x < xmax + pad ] < xmin - pad ]
		profile_fit = Spectrum( igauss(x, *coeff) )
		profile_fit.wave = x

		return profile_fit

	except OptionsError as err:
		print(' --> OptionsError:', err)
		raise ProfileError('Failed keyword assignment in Fit().')

	except DataTypeError as err:
		print(' --> DataTypeError:', err)
		raise ProfileError('DataTypeError in Fit().')

    
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
            'kind' : 'cubic', # given to scipy...interp1d for continuum
        })
        
        kind = options('kind')
        
    except OptionsError as err:
        print(' --> OptionsError:', err)
        raise ProfileError('Unrecognized option given to Extract()!')
        
    print(' Please select four points identifying the spectral line.')

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
    # keep wavelengths here for later
    xx = xc.copy()
    # but not the line
    yc = yc[np.where(np.logical_or(xc < wave[1], xc > wave[2]))]
    xc = xc[np.where(np.logical_or(xc < wave[1], xc > wave[2]))]
    
    # now use interpolation to cross the gap
    interp    = interp1d(xc, yc, kind = kind)
    continuum = Spectrum(interp(xx), xx)

    # display visual aid
    plt.plot(xx, interp(xx), 'r--')
    plt.fill_between(xl, yl, interp(xl), color = 'blue', alpha = 0.5)
    plt.draw()

    return line, continuum