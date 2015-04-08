#
# script for creating plot for NSF Proposal
#
def NSFplot1( *vargs, **kwargs ):

	import matplotlib, pyfits
	from   matplotlib  import pyplot as plot
	from   Pylib.astro import simbad, fits, cal as calibrate
	from   Pylib.astrolibpy.astrolib.helcorr import helcorr

	option = {}
	option['verbose']  = True
	option['barycorr'] = True

	for keyword in kwargs:

		if keyword not in option:
			raise ValueError("[%s] is not an available keyword argument!" % keyword)

		elif type( kwargs[keyword] ) is not type( option[keyword] ):
			raise TypeError("Keyword argument <%s> expects %s" % (keyword, type(option[keyword])))

		else: option[keyword] = kwargs[keyword]

	if len(vargs) > 0:

		if type(vargs[0]) is not list:
			raise TypeError("First variable argument must be a list of numpy arrays!")
		
		else: corrected = vargs[0]

	else:
		
		objects   = fits.getdata("ism/elodie/ten/raw",        verbose=option['verbose'])
		calibs    = fits.getdata("ism/elodie/allsky/regulus", verbose=option['verbose'])
		corrected = calibrate.telluric(objects, calibs,       verbose=option['verbose'])

	paths = fits.find("ism/elodie/ten/raw")

	marked = [3, 14, 20, 21, 22, 24]

	graph = [graph     for i, graph     in enumerate(corrected) if i in marked]
	paths = [this_path for i, this_path in enumerate(paths)     if i in marked]

	if option['barycorr']:

		if option['verbose']: print("\nPerforming Barycentric corrections...")

		# spatial coordinates of L'Observatoire de Haute-Provence (site of Elodie Spectrograph)
		latitude, longitude, altitude = 43.930833, 5.713333, 650.0 # deg, deg, meters
	
		for i, this_path in enumerate(paths):

			hdulist   = pyfits.open(this_path)[0]
			[ra, dec] = simbad.query( hdulist.header['OBJECT'], 'position')
			barycorr  = helcorr(latitude, longitude, altitude, ra*(24/360), dec, hdulist.header['JDB1'] )[0]
		
			graph[i][0] += graph[i][0] * barycorr / 299792.458	

	if option['verbose']: print("\nBuilding plot ...\n")

	matplotlib.rcParams['figure.facecolor'] = 'w'

	plot.ion()

	#fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plot.subplots(6, sharex=True, sharey=True)
	fig, (ax1, ax2, ax3, ax4) = plot.subplots(4, sharex=True, sharey=True)

	fig.set_size_inches(8, 4)

	ax1.plot(graph[0][0], graph[0][1], 'k-')
	ax2.plot(graph[2][0], graph[2][1], 'k-')
	ax3.plot(graph[4][0], graph[4][1], 'k-')
	ax4.plot(graph[5][0], graph[5][1], 'k-')
	#ax5.plot(graph[4][0], graph[4][1], 'k-')
	#ax6.plot(graph[5][0], graph[5][1], 'k-')

	plot.xlim(5885, 5900)
	plot.ylim(0.41, 1.17)

	ax1.tick_params(axis='both', which='major', labelsize=8)
	ax2.tick_params(axis='both', which='major', labelsize=8)
	ax3.tick_params(axis='both', which='major', labelsize=8)
	plot.tick_params(axis='both', which='major', labelsize=8)

	plot.rc('text', usetex=True)
	plot.rc('font', family='serif')

	fig.text(0.5, 0.01, "Wavelength (\AA)",   va='center', ha='center', size=12)
	fig.text(0.07, 0.5, "Relative Intensity", va='center', ha='center', size=12, rotation='vertical')

	path = fits.find("ism/elodie/ten/raw")

	path = [this_path for i, this_path in enumerate(path) if i in marked]

	identifier = []

	for this_path in path:

		hdulist = pyfits.open(this_path)

		identifier.append( hdulist[0].header['OBJECT'] )

	fig.text(0.22, 0.74, identifier[0], va='center', ha='center', size=10)
	#fig.text(0.22, 0.66, identifier[1], va='center', ha='center', size=10)
	fig.text(0.22, 0.54, identifier[2], va='center', ha='center', size=10)
	#fig.text(0.22, 0.40, identifier[3], va='center', ha='center', size=10)
	fig.text(0.22, 0.34, identifier[4], va='center', ha='center', size=10)
	fig.text(0.22, 0.139, identifier[5], va='center', ha='center', size=10)

	#fig.tight_layout()
	fig.subplots_adjust(hspace=0)

	plot.show()

def NSFplot2(*vargs, **kwargs):

	import numpy, matplotlib, pyfits
	from matplotlib import pyplot as plot
	from Pylib.astro import fits, simbad
	from Pylib.etc import display

	plot.ion()

	matplotlib.rcParams['figure.facecolor'] = 'w'

	option = {}
	option['verbose'] = True
	
	file_paths = fits.find("ism/elodie/ten/raw")

	if len(vargs) == 0:
	
		ra  = []
		dec = []
		dis = []

		if option['verbose']:
			print('Fetching positions & distances from SIMBAD ... ')

		for i, this_path in enumerate(file_paths):

			identifier = pyfits.open(this_path)[0].header['OBJECT']
			position   = simbad.query(identifier, 'position')
			distance   = simbad.query(identifier, 'distance')

			ra.append(  position[0] )
			dec.append( position[1] )
			dis.append( distance    )

			if option['verbose']: display.progress( (i+1)/len(file_paths) )
	
		ra  = numpy.array(ra)
		dec = numpy.array(dec)
		dis = numpy.array(dec)

	elif len(vargs) != 3:
		raise ValueError("must provide exactly three arguments if any!")

	else:
		ra = vargs[0]
		dec = vargs[1]
		dis = vargs[2]

	area = dis ** 2
	area_max = area.max( )

	area = 40.0 - ( 38.0 * area / area.max() + 1.0 ) # between 1 and 3 reversed

	plot.rc('text', usetex=True)
	plot.rc('font', family='serif')
	plot.figure(1, figsize=(7, 7))

	plot.clf()

	plot.scatter(ra, dec, s=area, c='k', alpha=1.0, label='data')
	
	plot.xlabel("Right Ascension", fontsize=14, labelpad=10)
	plot.ylabel("Declination",     fontsize=14, labelpad=10)

	# these are the index values for the spectra with ISM Na absorption
	marked = [2, 3, 12, 13, 14, 15, 17, 20, 21, 22, 23, 24, 26, 30, 32, 38, 39, 45]	

	target_ra  = [position for index, position in enumerate(ra)  if index in marked]
	target_dec = [position for index, position in enumerate(dec) if index in marked]
	target_area = [distance for index, distance in enumerate(area) if index in marked]
	
	target_ra = numpy.array(target_ra)
	target_dec = numpy.array(target_dec)
	target_area = 10.0 * numpy.array(target_area)

	plot.scatter(target_ra, target_dec, s=100, facecolors='none', 
	edgecolors='g', alpha=1.0, marker='D', label='detection')

	# these are the maybes
	marked = [1,6,7,8,10,11,16,19,25,27,29,31,33,36,37,40,46,47,50]

	target_ra   = [position for index, position in enumerate(ra)   if index in marked]
	target_dec  = [position for index, position in enumerate(dec)  if index in marked]
	target_area = [distance for index, distance in enumerate(area) if index in marked]
	
	target_ra   = numpy.array(target_ra)
	target_dec  = numpy.array(target_dec)
	target_area = 10.0 * numpy.array(target_area)

	plot.scatter(target_ra, target_dec, s=100, facecolors='none', 
	edgecolors='b', alpha=1.0, marker='s', label='potential')

	delcyg = simbad.query('del cyg', 'position')

	plot.scatter( delcyg[0], delcyg[1], s=120, facecolors='none',
	edgecolors='r', alpha=1.0, marker='o', label='$\delta$ Cygnus')

	plot.ylim(35, 56)

	plot.legend(loc='upper right', fontsize=11.5, ncol=4, markerscale=0.75, scatterpoints=1)

	distance = simbad.query('del cyg', 'distance') ** 2

	area = 40.0 - ( 38.0 * distance / area_max + 1.0 )

	plot.scatter( delcyg[0], delcyg[1], s=area, facecolors='r',
	edgecolors='r', alpha=1.0, marker='o')

	plot.grid(True)
		
	plot.show()

	if len(vargs) == 0:
		return ra, dec, dis
