# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# AstroPython/Astro/Observatory.py 
"""
Classes for defining observatory parameters (similar to the IRAF task).
"""

class Observatory:
	"""
	The Abstract base class for Observatory types.
	"""
	def __init__(self):
		raise TypeError('The Observatory base class should not be '
		'instantiated on its own.')

class OHP(Observatory):
	"""
	The Observatoire de Haute-Provence, France.
	"""
	def __init__(self):
		self.name      = 'Observatoire de Haute-Provence'
		self.latitude  = 43.9308334 # degrees N
		self.longitude = 356.28667  # degrees W
		self.altitude  = 650        # meters
