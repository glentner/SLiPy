# Copyright (c) Geoffrey Lentner 2015. All Right Reserved.
# See LICENSE (GPLv3)
# slipy/Algorithms/Functions.py
"""
Analytical functions used in SLiPy.
"""

import numpy as np
from scipy.special import wofz as w # Faddeeva function

# No FunctionsError implimented yet

def Gaussian(x, *params):
    """
    The Gaussian function. The function template has `*params` to work with
    scipy...curve_fit; the user *must* specify `A, mu, sigma = params`.
    """
    A, mu, sigma = params
    return A * np.exp( -0.5 * (x - mu)**2 / sigma**2 )

def NormalizedGaussian(x, *params):
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
    The Lorentzian function. Typically used to describe an absorption line profile.
    """
    x0, gamma = params
    return 1 / ((2 * (x - x0) / gamma)**2 + 1)

def NormalizedLorentzian(x, *params):
    """
    The Lorentzian function, normalized.
    """
    x0, gamma = params
    return 2 * Lorentzian(x, *params) / (np.pi * gamma)

def NormalizedVoigt(x, *params):
    """
    The Voigt function. The result of the convolution of a Lorentzian with one or more
    Gaussians. Often used to describe an intrinsically Lorentzian absorption profile blurred
    by Gaussian broadening from thermal motions and turbulance and another Gaussian instrument
    profile. The result is also the real part of the Faddevva function (the complex probability
    function), implemented by Scipy as wofz in scipy.special

    This returns a normalized profile. The name `NormalizedVoigt` is to keep the convention
    for the rest of this module. The amplitude controlling Voigt profile `Voigt` evaluates
    this function first for `x - x0 = 0` to empirically `un-normalize` it before scaling by the
    requested amplitude

    Reference:
        Olivero, J.J; R.L. Longbothum. ``Empirical fits to the Voigt line width:
        A brief review''. Journal of Quantitative Spectroscopy and Radiative Transfer. Vol 17.
        Issue 2. pages 223-236. Feb 1977

    Parameters: `x0, sigma, gamma = params`

        x0    -> center of the profile
        sigma -> the Gaussian width
        gamma -> the Lorentzian width
    """
    x0, sigma, gamma = params
    return w( ((x-x0) + 1j*gamma) / (sigma * np.sqrt(np.pi)) ).real / (
        sigma * np.sqrt(2 * np.pi) )

def Voigt(x, *params):
    """
    The Voigt line profile. See ..Algorithms.Function.NormalizedVoigt for more information.
    This function returns an amplitude controlled Voigt profile.

    Parameters: `A, x0, sigma, gamma = params`

        A     -> amplitude of the profile
        x0    -> center of the profile
        sigma -> the Gaussian width
        gamma -> the Lorentzian width
    """
    return params[0] * NormalizedVoigt(x, *params[1:]) / NormalizedVoigt(0, 0, *params[2:])


def InvertedLorentzian(x, *params):
    """
    An inverted Lorentzian (i.e, A - Lorentzian()).
    """
    A, x0, gamma = params

    return A - Lorentzian(x, x0, gamma)
