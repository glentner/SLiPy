# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# SLiPy/__init__.py
"""
SLiPy
===========
A Computational Astronomy Package for Python

This Python package contains my code for performing computational work in
astronomy research.
"""

# TkAgg must be used first for Mac user's
# import matplotlib as mpl
# mpl.use('TkAgg')

# exposed modules
from .SLiPy import Fits, Simbad, Correlate, Telluric, Velocity, \
		Observatory, Montage, Plot, DataType, Measure

# from .SLiPy.Plot     import SPlot, Iterate
# from .SLiPy.DataType import Spectrum
# from .SLiPy.Measure  import Deblend

# from .Data import Elodie

from .Framework import Display
