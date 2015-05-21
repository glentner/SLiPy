# Copyright (c) Geoffrey Lentner 2015. All Right Reserved.
# See LICENSE (GPLv3)
# slipy/Algorithms/Functions.py
"""
Analytical functions used in SLiPy.
"""

import numpy as np

# No FunctionsError implimented yet

def Gaussian(x, *params):
    """
    The Gaussian function. The function template has `*params` to work with
    scipy...curve_fit; the user *must* specify `A, mu, sigma = params`.
    """
    A, mu, sigma = params
    return A * np.exp( -0.5 * (x - mu)**2 / sigma**2 )

def NormalizedGaussian(x *params):
    """
    The normalized Gaussian function.
    """
    mu, sigma = params
    return np.exp( -0.5 * (x - mu)**2 / sigma**2 ) / (sigma *
        np.sqrt(2 * np.pi))

def InvertedGaussian(x, *params):
    """
    An inverted Gaussian function (i.e., 1 - Guassian()). The function
    template has `*params` to work with scipy...curve_fit; the user *must*
    specify `A, mu, sigma = params`.
    """
    return 1 - Gaussian(x, *params)

def Lorentzian(x, *params):
    """
    The Lorentzian function.
    """
    x0, gamma = params
    return (0.5 * gamma / np.pi) / ( (x - x0)**2 + (0.5 * gamma)**2 )

def InvertedLorentzian(x, *params):
    """
    An inverted Lorentzian (i.e, A - Lorentzian()).
    """
    A, x0, gamma = params

    return A - Lorentzian(x, x0, gamma)
