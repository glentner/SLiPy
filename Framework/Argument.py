# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/Framework/Arguments.py
"""
Module contains `Argument` class for handling conversions and type
checking for function/class keyword argument options.
"""

from .. import SlipyError

class ArgumentError(SlipyError):
	"""
	Exception specific to Argument module.
	"""
	pass

class Argument:
	"""
	`Argument` object has type and value management.
	"""
	def __init__(self, value, name = 'unspecified', **kwargs ):
		"""
		Build Argument `value`, `type`, and set `options`.
		"""
		options = {
				'lock' : False # disallow all type conversions
			}
		for arg in kwargs:
			if arg not in options:
				raise ArgumentError('`{}` is not an option for Argument.'
					.format(arg))
			if type(kwarg[arg]) is not type(options[arg]):
				raise ArgumentError('Option `{}` expects {}, given {}'
					.format(arg, type(options[arg]), type(kwargs[arg])) )
			# accept reassignment
			options[arg] = kwargs[arg]

		self.lock  = options['lock']
		self.T     = type(value)
		self.value = value
		self.name  = str(name)
		if type(name) is not str:
			raise ArgumentError('Argument expects type str for name.')

	def __call__(self, value):

		if self.lock and type(value) is not self.T:
			raise ArgumentError('Argument `{}` is locked as {}, '
				'but new value has {}.'
				.format(self.name, self.T, type(value)) )

		if self.T is bool:
			# special rules for `bool` conversions

			if type(value) is str:
				if value != "True" and value != "False":
					raise ArgumentError('Invalid conversion from {} '
						'for Argument `{}`.'.format(self.T, self.name))
				self.value = True if value == 'True' else False

			elif type(value) is int:
				if value != 0 and value != 1:
					raise ArgumentError('Invalid conversion from {} '
						'for Argument `{}`.'.format(self.T, self.name) )
				else:
					self.value = True if value else False

			elif type(value) is not bool:
				raise Error('Invalid conversion from {} '
					'for Argument `{}`.'.format(self.T, self.name))

			else: self.value = value

		else:

			try:
				self.value = self.T(value)

			except ValueError as error:
				raise ArgumentError('Cannot convert {} to {} for `{}`'
					' Argument.'.format(self.T, type(value), self.name))
