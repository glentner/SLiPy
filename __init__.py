# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/__init__.py
"""
AstroPython
===========
A Computational Astronomy Package for Python

This Python package contains my code for performing computational work in
astronomy research.
"""

# exposed modules
from .AstroPython import Fits, Simbad, Correlate, Telluric, Velocity, \
		Observatory, Montage

from .AstroPython.Plot     import SPlot, Iterate
from .AstroPython.DataType import Spectrum

from .Data import Elodie

from .Framework import Display
