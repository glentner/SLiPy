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

from ..Algorithms.Functions import Gaussian, InvertedLorentzian
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

def Fit(splot, function = InvertedLorentzian, params = None):
    """
    Given `splot` of type Plot.SPlot, the user selects two points on the 
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

def Remove(spectrum, line, continuum, **kwargs):
    """
    Given a continuum normalized `spectrum`, that contains a `line` profile
    that has been modelled with it's surrounding `continuum`, remove the
    line from the spectrum by division.
    """