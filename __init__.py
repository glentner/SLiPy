# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# SLiPy/__init__.py
"""
SLiPy - A Spectroscopy and astrophysics Library for Python 3

This Python package is an expanding code base for doing computational
astronomy, particularly spectroscopy. It contains both a *Spectrum* class
for handling spectra as objects (with +, -, \*, /, etc... operations defined)
and a growing suite of analysis tools.
"""

# base exception class for whole project, module exception classes will
# be derived from here.
class SlipyError(Exception):
    pass

# exposed modules
from .SLiPy import Fits, Simbad, Correlate, Telluric, Velocity, \
		Observatory, Montage, Plot, DataType, Measure, Profile

# elevate `Spectrum` to the package level
from .SLiPy.DataType import Spectrum