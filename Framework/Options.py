# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/Framework/Options.py
"""
Class object for handling kwargs in classes and functions.
"""

from .. import SlipyError
from .Argument import Argument as Arg, ArgumentError

class OptionsError(SlipyError):
	"""
	Exception specific to Options module.
	"""
	pass

class Options:
	"""
	Class object for handling kwargs in classes and functions.
	"""
	def __init__(self, kwargs, options ):
		"""
		Check types and build options dictionary
		"""
		try:
			# initial assignment
			self.options = {
					name : Arg(value, name)
					for name, value in options.items()
				}
			# attempted reassignment
			for key, value in kwargs.items():
				self.options[key](value)

		except AttributeError as err:
			raise OptionsError(
				'Options object expects dictionary types.')

		except KeyError as key:
			raise OptionsError(
				'{} was not a recognized option.'.format(key))

		except ArgumentError as err:
			print('\n --> ArgumentError:', err.msg )
			raise OptionsError('Failed assignment.')

	def __call__(self, option):
		"""
		Retrieve value of `option`.
		"""
		try:
			return self.options[option].value
		except KeyError as key:
			raise OptionsError('{} was not recognized.'.format(key))

	def items(self):
		"""
		Access options.items() values.
		"""
		return { k:v.value for k,v in self.options.items() }.items()
