# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Python/__init__.py
"""
Python Library (of Geoffrey Lentner)
"""

class BaseError(Exception):
	"""
	Base Exception class for the Python Library 
	"""
	def __init__(self, msg = 'undefined', arg = None):
		"""
		Construct the error with a message and an optional argument
		"""
		self.msg = msg
		self.arg = arg

# all subpackages
__all__ = ['Astro','astrolibpy','General']
