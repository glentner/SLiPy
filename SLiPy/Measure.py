# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# slipy/SLiPy/Measure.py

'''
Measure - Collection of measurement/calculation tools for spectra, etc.
'''
# import as TkAgg for widget tools

from matplotlib.widgets import  RectangleSelector
from pylab import *

from .Plot import SPlot, PlotError
from .DataType import Spectrum, DataError

class MeasureError(Exception):
    '''
    Exception specific to the Measure module.
    '''
    pass

#
# The below widgets were adapted from the example usage at
# matplotlib.org/api/widgets_api.html
#

# global variables accessable to the `selector` functions
regions = [ ]
plot = None

def OnSelect_Deblend(eclick, erelease):
    '''Selector function specific to the Deblend() function.'''

    # region selected
    x1, x2 = eclick.xdata, erelease.xdata

    # make sure we are in the correct order
    if x2 < x1:
        x1, x2 = x2, x1

    # grab `regions` list
    global regions
    global plot

    # x1 *= plot.wave[0].unit
    # x2 *= plot.wave[0].unit

    print(x1, x2)

    # create new Spectrum of region
    # data = plot.data[0].copy()
    # wave = plot.wave[0].copy()
    # data = data[ wave[wave < x2] > x1 ]
    # wave = wave[ wave[wave < x2] > x1 ]
    #
    # this_region = Spectrum(data, wave)

    # append selected region to list
    # regions.append( SPlot(this_region) )

    # add it to the current `SPlot`
    # plot.overlay( regions[-1] )

    # make it a red dashed line
    # plot.marker[-1] = 'r-'

    # push updates
    # plot.refresh()

def toggle_selector(event):
    '''activate/deactivate selection tool'''
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' Selector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' Selector activated.')
        toggle_selector.RS.set_active(True)

def Deblend(spectrum, **kwargs):
    '''
    Remove feature from `spectrum` via subtraction. The user manually selects
    regions in the `spectrum` and a `line` is fit to and subtracted. The
    original `spectrum` object is modified.

    If `spectrum` is actually a Spectrum object (DataType.Spectrum), a `SPlot`
    figure is created for it (without a label). Alternatively, a SPlot object
    can be directly given.
    '''
    if type(spectrum) is Spectrum:
        # create a SPlot object out of the spectrum
        spectrum = SPlot(spectrum)
    elif type(spectrum) is not SPlot:
        raise MeasureError('From Deblend(), argument was be either of type '
        'Spectrum or SPlot!')

    # make spectrum window active
    spectrum.draw()

    # pass off spectrum to global plot object
    global plot
    plot = spectrum


    toggle_selector.RS = RectangleSelector(plot.ax, OnSelect_Deblend,
    	minspanx=1, minspany=1)

    connect('key_press_event', toggle_selector)

    print('Select regions for continuum fitting: <press enter when done>')
    if input(): return plot
