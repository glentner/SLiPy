# Copyright (c) Geoffrey Lentner 2015. All Right Reserved.
# See LICENSE (GPLv3)
# slipy/Algorithms/KernelFit.py
"""
Non-parametric Kernel Regression.
"""

import numpy as np

from .. import SlipyError
from ..Framework.Options import Options, OptionsError

from .Functions import Gaussian

class KernelFitError(SlipyError):
    """
    Exception specific to the KernelFit module
    """

class KernelFit1D():
    """
    One dimensional kernel regression.
    """
    def __init__(self, x, y, kernel = Gaussian, **kwargs):
        """
        Keep the input `x` and `y` arrays. Define parameters.
        The default kernel function is the Guassian. Alternatives should
        have the same signature as slipy.Algorithms.Functions.Gaussian
        """
        if not hasattr(x, '__iter__') or not hasattr(y, '__iter__'):
            raise KernelFitError('Both `x` and `y` are expected to '
            'be array-like!')
        
        if len(np.shape(x)) != 1 or len(np.shape(y)) != 1:
            raise KernelFitError('Both `x` and `y` are expected to '
            'be one-dimensional objects!')
        
        if len(x) != len(y):
            raise KernelFitError('Both `x` and `y` should be of the '
            'same length!')
            
        if not hasattr(x, 'copy') or not hasattr(y, 'copy'):
            raise KernelFitError('Both `x` and `y` should have a '
            'copy() method implimented!')
        
        if type(x) is not type(y):
            raise KernelFitError('Both `x` and `y` are expected to '
            'be of the same type!')
        
        if hasattr(x, 'unit') and not hasattr(y, 'unit'):
            raise KernelFitError('The `x` array given has units but the `y` '
            'array doesn`t!, maybe give `y` u.dimensionless_unscaled?')
        
        if hasattr(y, 'unit') and not hasattr(x, 'unit'):
            raise KernelFitError('The `y` array given has units but the `x` '
            'array doesn`t!, maybe give `x` u.dimensionless_unscaled?')
            
        # the default `bandwidth` is 1/10 the domain of `x`, this
        # is likely to be a bad choice!
        bandwidth = 0.1 * (x.max() - x.min())
        
        try:
            
            options = Options( kwargs, {
                
                'bandwidth' : bandwidth # for the kernel function
            })
            
            self.kernel    = kernel
            self.bandwidth = options('bandwidth') # not necessarily with units
            
        except OptionsError as err:
            print(' --> OptionsError:', err)
            raise KernelFitError('Unrecognized option for KernelFit1D()!')
        
        self.x = x.copy()
        self.y = y.copy()
    
    def mean(self, x):
        """
        Solve for smooth profile through the data on the new `x` array.
        This is essentially a weighted mean.
        """
        
        if not hasattr(x, '__iter__'):
            raise KernelFitError('`x` is expected to be array-like!')
        
        if not hasattr(x, 'copy'):
            raise KernelFitError('`x` is expected to have a copy() method!')
        
        if hasattr(x, 'unit') and not hasattr(self.x, 'unit'):
            raise KernelFitError('The provided array has units but the '
            'original domain does not!')
        
        if hasattr(self.x, 'unit') and not hasattr(x, 'unit'):
            raise KernelFitError('The provided array does not have units '
            'but the original domain did!')
        
        # copy `x` such that the return type is the same as the input type
        y = x.copy()
        
        # fix units though
        if hasattr(x, 'unit'):
            y = y.value * self.y.unit
        
        for a, point in enumerate(x):
            
            weights = Gaussian(point, 1, self.x, self.bandwidth)
            
            # weights should be dimensionless
            if hasattr(x, 'unit'):
                weights = weights.decompose()
            
            y[a] = np.sum(weights * self.y) / np.sum(weights)
        
        return y
        