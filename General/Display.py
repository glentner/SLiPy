# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv2)
# Display/General/Display.py
"""
Display.py - Python module for displaying content to the screen.
"""

import os, sys, math
from time import time
from datetime import datetime, timedelta

from Python import BaseError
from Python.General.Options import *

class DisplayError(BaseError):
	"""
	Module specific exception.
	"""
	pass

class Monitor:
	"""
	Class for displaying a progress bar during iterative tasks.
	"""
	def __init__(self, **kwargs ):
		try:
			# available keyword options 
			self.options = Options( kwargs,
				{
					'width'    : 45    , # number of characters wide
					'numbers'  : True  , # display numberical percent
					'template' : '[=>]', # template for progress bars
					'freq'     : 0.25  , # refresh rate
					'ETC'      : False , # display estimated time of completion
					'inline'   : True    # vanish after completion
				})

			# give assignments
			self.width   = self.options('width')
			self.numbers = self.options('numbers')
			self.freq    = self.options('freq')
			self.ETC     = self.options('ETC')
			self.inline  = self.options('inline')
			self.left, self.char, self.tip, self.right = self.options('template')

			# start clocks
			self.start = time()
			self.last  = time()

		except OptionsError as err:
			print( '\n --> OptionsError:', err.msg )
			raise DisplayError('Failed to initialize Monitor.')

		except ValueError as err:
			raise DisplayError(
				'`template` option requires exactly 4 characters.')

	def __EstimatedCompletionTime(self):
		"""
		Estimated time of completion, based off percent complete and
		current total elapsed time.
		"""
		if self.percent > 0:
			elapsed   = time() - self.start
			remaining = elapsed * (1 / self.percent - 1)
			etc = datetime.today() + timedelta(seconds=remaining) 
			return etc.strftime(' ETC: %Y-%m-%d @ %H:%M ') 

	def __build(self):
		"""
		Build the progress bar
		"""
		bars   = self.char * math.floor( self.percent * self.width )
		empty  = ' ' * (self.width - len(bars) - 1)
		display = self.left + bars + self.tip + empty + self.right 

		if self.numbers:
			display += '{:>8.2f} % '.format( self.percent * 100 )

		if self.ETC:
			display += self.__EstimatedCompletionTime()
		
		sys.stdout.write('\r \033[K \r {}'.format(display))
		sys.stdout.flush()

	def progress(self, i, imax):
		"""
		Request a progress bar.
		"""
		if time() - self.last > self.freq:
			# refresh rate surpasses,
			# update time of last call and percent complete
			self.last    = time()
			self.percent = float(i) / float(imax) 
			# display progress bar
			self.__build()

	def complete(self):
		"""
		Call to finalize the progress bar.
		"""
		if self.inline:
			sys.stdout.write('\r\033[K\r')
			sys.stdout.flush()

		else:
			self.percent = 1
			self.numbers = False
			self.ETC     = False
			self.__build()
			sys.stdout.write(' complete\n')
			sys.stdout.flush()

	def elapsed(self):
		"""
		Display total time elapsed since instantiation.
		"""
		total  = time() - self.start 
		abrv  = [ 'd',        'h',      'm',    's']
		unit  = { 'd': 86400, 'h':3600, 'm':60 }
		count = { 'd':0,      'h':0,    'm':0,  's':0 }

		for item in abrv:
			if item in unit:
				while total > unit[item]:
					total -= unit[item]
					count[item] += 1
			else: count[item] = math.floor(total)

		total = [	
				'{} {}'.format( v, u )
				for u, v in zip( abrv, [ count[v] for v in abrv ] ) 
				if count[u]
			]
		total = ' Time Elapsed: ' + ' '.join(total)
		total = ' ' + '-' * (len(total)+5) + '\n ' + total

		sys.stdout.write(total)
		sys.stdout.flush()
