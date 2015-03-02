# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/__init__.py
"""
 - Computational Astronomy Package for Python -

This Python package contains my code for performing computational work in 
astronomy research. I've originally named this package simply `Python` because
there seems to be an over abundance of acronyms in the current state of 
astronomy computing with Python that I didn't want to pollute it any further.
As such, I've implemented this project with `relative` imports within the 
package source such that you can rename it whatever you please (i.e.,
git clone http://github.com/glentner/Python; scp -r Python Susy)

"""

# exposed modules
from .Astro import Fits, Simbad, Plot, DataType
