#
#
#

def setup():

	from matplotlib import pyplot as plot
	import matplotlib

	matplotlib.rcParams['figure.facecolor'] = 'w'
	plot.rc('text', usetex=True)
	plot.rc('font', family='serif')
	plot.ion()

def splot( data, *vargs, **kwargs ):

	from   matplotlib import pyplot as plot
	import matplotlib
	import numpy

	matplotlib.rcParams['figure.facecolor'] = 'w'

	option = {}
	option['fancy']   = False
	option['iterate'] = False
	option['marker']  = 'k-'

	for argument in vargs:

		if argument not in option:
			raise ValueError( "[%s] is not a recognized option!" % argument )

		elif type( option[argument] ) is not bool:
			raise TypeError( "[%s] is not available for assignment this way!" % argument )

		else: option[argument] = True

	for keyword in kwargs:

		if keyword not in option:
			raise ValueError( "[%s] is not a recognized option!" % argument)

		elif type( kwargs[keyword] ) is not type( option[keyword] ):
			raise TypeError( "[%s] must have be of %s" % (keyword, type(option[keyword])) )

		else: option[keyword] = kwargs[keyword]
	
	if option['fancy']:
		plot.rc('text', usetex=True)
		plot.rc('font', family='serif')
	
	plot.ion()
			
	mark = [];

	if type( data ) is list:
		
		for a, graph in enumerate(data):
			
			if option['iterate']: 
				
				plot.clf()
				plot.plot(graph[0], graph[1], option['marker'])
				plot.title("Graph %s out of %s" % (a+1, len(data)))
				plot.xlabel("Wavelength (Angstroms)")
				plot.ylabel("Normalized Intensity")
				plot.xlim(graph[0][0], graph[0][-1])

				response = input("Press <return> to continue (press any key to save index):")

				if response != "": mark.append(a)

			else: plot.plot(graph[0], graph[1])

	else: plot.plot(data[0], data[1], option['marker'])

	plot.xlabel("Wavelength (Angstroms)")
	plot.ylabel("Normalized Intensity")

	return mark
