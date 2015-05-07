# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2).
# slipy/Framework/Interface.py 
"""
Interface.py

Command Line Interface tools.
"""

from .Options import Options, OptionsError

class CommandError(Exception):
	"""
	Exception for Interface module.
	"""
	pass

def Parse( clargs, **kwargs ):
	"""
	Parse command line arguments, `clargs` (i.e., sys.argv).
	"""
	if type(clargs) is not list:
		raise CommandError('Parse function expects a list for `clargs`.')
	try:
		options = Options( kwargs,
			{
				'exe' : True # search first argument for '@func' pattern
			})

		if options('exe'):
			function = clargs[0].split('@')
			if len(function) != 2 or function[0]:
				raise CommandError('Incorrect formatting in function call.')
			function = function[1]
			del(clargs[0])

		# args should not have an assignment
		args = [ x for x in clargs if '=' not in x ]

		# remaining clargs should be kwargs
		kwargs = {
				key : value for key, value in [
					arg.split('=') for arg in set(clargs) - set(args)
					]
			}

		if options('exe'):
			return function, args, kwargs

		else:
			return args, kwargs

	except OptionsError as err:
		print('\n --> OptionsError:', err)
		raise CommandError('from Parse')

	except ValueError as key:
		raise CommandError('Incorrect formatting of keyword arguments.')
