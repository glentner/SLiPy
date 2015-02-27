#
# File:        cal.py
# Description: Module contains function definitions for calibration tasks on spectra
#
# Copyright (C) Geoffrey Lentner 2014
# See LICENCE under main directory
#
# Contact: Geoffrey Lentner
#          Graduate Student / Researcher
#          LL15 Natural Science Building
#          Department of Physics & Astronomy
#          University of Louisville
#          Louisville, KY 40291 USA
#          geoffrey.lentner@louisville.edu
#
# Update: 10.6.2014 @ 21:15 EST
#
 
# Routine for x-correlation spectra  with specified lag
def correlate ( array_1, array_2, **kwargs ):

	import numpy
	
	from Pylib.etc import display

	# define options for function
	option            = {}
	option['verbose'] = False
	option['lag']     = 25

	for item in kwargs:
		if item not in option:
			raise ValueError( "[%s] is not a recognized option!" )

	if 'verbose' in kwargs:

		option['verbose'] = kwargs['verbose']

		if type( option['verbose'] ) is not bool:
			raise TypeError( "Option [verbose] must be of type bool!" )
	
	if 'lag' in kwargs:

		option['lag'] = kwargs['lag']

		if type( option['lag'] ) is not int:
			raise TypeError( "Option [lag] must be of type int!" )

	if type( array_1 ) is not numpy.ndarray or type( array_2 ) is not numpy.ndarray:
		raise TypeError( "Input arrays should be of type numpy.ndarray!" )

	if len( numpy.shape( array_1 ) ) != 1 or len( numpy.shape( array_2 ) ) != 1:
		raise ValueError( "Input arrays should be 1D numpy.ndarray's!" )

	if len( array_1 ) != len( array_2 ):
		raise ValueError( "Input arrays should be of the same length!" )

	if option['lag'] > len( array_1 ) - 1:
		raise ValueError( "Option [lag] must be at least one less than length of spectra!" )

	if option['verbose']: print("\nCorrelating arrays ... ")

	size = len( array_1 )
	
	rms = []

	for i, shift in enumerate( range( -option['lag'], option['lag'] + 1 ) ):

		if    shift < 0: difference = array_1[-shift:] - array_2[:shift]
		elif  shift > 0: difference = array_1[:-shift] - array_2[shift:]
		else:            difference = array_1          - array_2

		rms.append( numpy.sqrt(( difference ** 2 ).sum( ) / size) )

		if option['verbose']: display.progress( i / (2*lag + 1) ) 

	return numpy.array(rms)


# Perform telluric correction on one or more spectra using one or more calibration spectra  #
def telluric ( OBJECT, CALIBRATION, **kwargs ):

	import os, numpy 
	
	from Pylib.astro import fits
	from Pylib.etc   import display
	from sys         import stdout; 
	
	write = stdout.write
    
	# define a default parameter dictionary              #
	option            = {}                               #
	option['verbose'] = True                             # provide progress statements
	option['xcorr']   = True                             # perform the cross correlations (horizontal)
	option['lag']     = 25                               # pixels to shift in xcorrelation
	option['amplify'] = True                             # perform the amplification fit (vertical)
	option['range']   = numpy.linspace(0.5, 2.0, 31)     # range of values to try
	option['trials']  = len( option['range'] )           # number of trials in 'range'

	for item in kwargs:
		if item not in option:
		    raise ValueError( "Argument [%s] not recognized!" % item )
	
	if 'xcorr' in kwargs:
        
		option['xcorr'] = kwargs['xcorr']

		if type( option['xcorr'] ) is not bool:
			raise TypeError ("Parameter [xcorr] must be of type bool")

	if 'lag' in kwargs:
    
		option['lag'] = kwargs['lag']
    
		if type( option['lag'] ) is not int:
			raise TypeError ( "Parameter [lag] must be of type int" )

	if 'amplify' in kwargs:

		option['amplify'] = kwargs['amplify']
        
		if type( option['amplify'] ) is not bool:
			raise TypeError ( "Parameter [amplify] must be of type bool" )

	if 'range' in kwargs:
        
		option['range']  = kwargs['range']
		option['trials'] = len( option['range'] )
        
		if type( option['range'] ) is not numpy.ndarray:
			raise TypeError ( "Parameter [range] must be of type numpy.ndarray" )

		elif numpy.shape( option['range'] )[0] != 1:
			raise ValueError ( "[range] must be a 1D array!" )

	if type( OBJECT ) is str: 
		OBJECT = fits.getdata( OBJECT, verbose=option['verbose'] )
	
	if type( CALIBRATION ) is str: 
		CALIBRATION = fits.getdata( CALIBRATION, verbose=option['verbose'] )
    
	if type( OBJECT ) is list:
		
		for i, item in enumerate( OBJECT ):
		
			if type( item ) is not numpy.ndarray:
				raise TypeError ( "Item %d in OBJECT list was not of type numpy.ndarray!" % i )

			if len( numpy.shape( item ) ) != 2:
				raise ValueError( "Item %d in OBJECT list must be 2D numpy.ndarray!" % i )

			if numpy.shape( item )[0] != 2:
				raise ValueError( "Item %d in OBJECT list must have exactly 2 rows!" % i )

	elif type( OBJECT ) is numpy.ndarray:

		if len( numpy.shape( OBJECT ) ) != 2:
			raise ValueError( "Object spectra must be 2D numpy.ndarray!" )

		if numpy.shape( OBJECT )[0] != 2:
			raise ValueError( "Object spectra must have exactly 2 rows!" )

	else: raise TypeError (" OBJECT spectra must be of type numpy.ndarray" )

	if type( CALIBRATION ) is list:
        
		for i, item in enumerate( CALIBRATION ):
            
			if type( item ) is not numpy.ndarray:
				raise TypeError( "Item %d in CALIBRATION list was not of type numpy.ndarray!" % i )
            
			if len( numpy.shape( item ) ) != 2:
				raise ValueError( "Item %d in CALIBRATION list must be 2D numpy.ndarray!" % i )

			if numpy.shape( item )[0] != 2:
				raise ValueError( "Item %d in CALIBRATION must have exactly 2 rows!" % i )

	elif type( CALIBRATION ) is numpy.ndarray:
		
		if len( numpy.shape( CALIBRATION ) ) != 2:
			raise ValueError( "Calibration spectrum must be 2D numpy.ndarray!" )

		if numpy.shape( CALIBRATION )[0] != 2:
			raise ValueError( "Calibration spectrum must have exactly 2 rows!" )
		
	else: raise TypeError( "CALIBRATION spectra must be of type numpy.ndarray!" )
		
	if type( OBJECT ) is list:

		comparison = numpy.shape( OBJECT[0] )[1]

		for spectra in OBJECT:
			if numpy.shape( spectra )[1] != comparison:
				raise ValueError( "Not all spectra provided are of equal length!" )

	else: comparison = numpy.shape( OBJECT )[1]

	if type( CALIBRATION ) is list:

		for spectra in CALIBRATION:
			if numpy.shape( spectra )[1] != comparison:
				raise ValueError( "Not all spectra provided are of equal length!" )

	else:
		if numpy.shape( CALIBRATION )[1] != comparison:
			raise ValueError( "Not all spectra provided are of equal length!" )

	if type( OBJECT ) is list: tasks = [ len( OBJECT ) ]
	else:                      tasks = [1]

	if type( CALIBRATION ) is list: tasks.append( len( CALIBRATION ) )
	else:                           tasks.append( [1] )

	if comparison * option['trials'] > 10**7:
		
		warning = "The matrices involved here have > 10^7 elements, proceed anyway? yes/[no]: "

		while True:

			feedback = input( warning )

			if   feedback == '' or feedback == 'n' or feedback == 'no':  return
			elif                   feedback == 'y' or feedback == 'yes': break
			else: print( "[%s] is not a recognized response!" % feedback )

	if not option['amplify'] and type( CALIBRATION ) is list:

		warning  = "Option [amplify] set to False but user provided CALIBRATION list. \n"
		warning += "There is no criteria for quantifying which calibration spectra is best.\n"
		
		raise ValueError( warning )

	if option['verbose']: 
		write( "\nPerforming telluric corrections on " )
		write( "%d object spectra with %d calibration spectra ... \n" % (tasks[0], tasks[1]) )	

	if type( OBJECT ) is list:

		for i, obj_spectra in enumerate( OBJECT ):

			if type( CALIBRATION ) is list:

				for j, cal_spectra in enumerate( CALIBRATION ):

					if option['xcorr']:
						
						shift = correlate(obj_spectra[1], cal_spectra[1], lag=option['lag'])
						shift = shift.argmin( ) - option['lag']

						if shift < 0:

							cal_matrix = numpy.tile( cal_spectra[1][:shift ], (option['trials'], 1) )
							obj_matrix = numpy.tile( obj_spectra[1][-shift:], (option['trials'], 1) )

						if shift > 0: 
							
							cal_matrix = numpy.tile( cal_spectra[1][shift: ], (option['trials'], 1) )
							obj_matrix = numpy.tile( obj_spectra[1][:-shift], (option['trials'], 1) )

						else:

							cal_matrix = numpy.tile( cal_spectra[1], (option['trials'], 1) )
							obj_matrix = numpy.tile( obj_spectra[1], (option['trials'], 1) )

						size       = numpy.shape( cal_matrix )[1]
						amplitude  = numpy.tile(option['range'], (size, 1)).T
						cal_matrix = ( 1 - ( 1 - cal_matrix ) * amplitude )
							
						difference  = obj_matrix - cal_matrix

						rms = numpy.sqrt( (difference ** 2).sum(axis=1) / size )
						
						if j == 0:

							best = rms.min( ); site = [rms.argmin( ), shift, j]

						elif best > rms.min( ): 
							
							best = rms.min( ); site = [rms.argmin( ), shift, j]
						
					else: raise ValueError( "Not coded yet!" ) 
		
				j = site[2]; shift = site[1]; site = site[0]

				if   shift < 0: OBJECT[i][1][-shift:] /= (1 - (1 - CALIBRATION[j][1][:shift]) * option['range'][site])
				elif shift > 0: OBJECT[i][1][:-shift] /= (1 - (1 - CALIBRATION[j][1][shift:]) * option['range'][site])
				else:           OBJECT[i][1]          /= (1 - (1 - CALIBRATION[j][1]        ) * option['range'][site])

			else: raise ValueError( "Not coded yet!" )
			
			if option['verbose']: display.progress( (i+1)/len(OBJECT) )

	else: raise ValueError( "Not coded yet!" )

	return OBJECT

