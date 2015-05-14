# Copyright (c) Geoffrey Lentner 2015. All Right Reserved.
# See LICENSE (GPLv3)
# slipy/Algorithms/Functions.py
"""
Analytical functions used in SLiPy.
"""

import numpy as np

# No FunctionsError implimented yet

def Gaussian(x, A = 1, mu = 0, sigma = 1):
    """
    The Gaussian function.
    """
    return A * np.exp( -0.5 * (x - mu)**2 / sigma**2 )
    
def InvertedGaussian(x, A = 1, mu = 0, sigma = 1):
    """
    An inverted Gaussian function (i.e., 1 - Guassian()).
    """
    return 1 - Gaussian(x, A, mu, sigma)