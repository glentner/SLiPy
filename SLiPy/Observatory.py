# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/SLiPy/Observatory.py
"""
Classes for defining observatory parameters (similar to the IRAF task).
"""

from astropy import units as u

# there is no need for an `ObservatoryError` class ...

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
        self.longitude = 356.28667 * u.degree # West
        self.latitude  = 43.9308334 * u.degree # North
        self.altitude  = 650 * u.meter
        self.timezone  = 1 * u.hourangle


#
# All of the below observatory parameters have been taken directly from IRAF!
#

class KPNO(Observatory):
	"""
	Kitt Peak National Observatory
	"""
	def __init__(self):
		self.name      = "Kitt Peak National Observatory"
		self.longitude = 111.6 * u.degree # West
		self.latitude  = 31.9633333333 * u.degree # North
		self.altitude  = 2120. * u.meter
		self.timezone  = 7 * u.hourangle

class WIYN(Observatory):
	"""
	WIYN Observatory
	"""
	def __init__(self):
		self.name      = "WIYN Observatory"
		self.longitude = 111.6 * u.degree # West
		self.latitude  = 31.9633333333 * u.degree # North
		self.altitude  = 2120. * u.meter
		self.timezone  = 7 * u.hourangle

class CTIO(Observatory):
	"""
	Cerro Tololo Interamerican Observatory
	"""
	def __init__(self):
		self.name      = "Cerro Tololo Interamerican Observatory"
		self.longitude = 70.815 * u.degree # West
		self.latitude  = -30.16527778 * u.degree # North
		self.altitude  = 2215. * u.meter
		self.timezone  = 4 * u.hourangle

class LASILLA(Observatory):
	"""
	European Southern Observatory: La Silla
	"""
	def __init__(self):
		self.name      = "European Southern Observatory: La Silla."
		self.longitude = 70.73 * u.degree # West
		self.latitude  = -28.7433333333 * u.degree # North
		self.altitude  = 2347 * u.meter
		self.timezone  = 4 * u.hourangle

class PARANAL(Observatory):
	"""
	European Southern Observatory: Paranal
	"""
	def __init__(self):
		self.name      = "European Southern Observatory: Paranal"
		self.longitude = 70.4033333333 * u.degree # West
		self.latitude  = -23.375 * u.degree # North
		self.altitude  = 2635 * u.meter
		self.timezone  = 4 * u.hourangle

class LICK(Observatory):
	"""
	Lick Observatory
	"""
	def __init__(self):
		self.name      = "Lick Observatory"
		self.longitude = 121.636666667 * u.degree # West
		self.latitude  = 37.3433333333 * u.degree # North
		self.altitude  = 1290 * u.meter
		self.timezone  = 8 * u.hourangle

# Observatory entry from a conversation with Craig Foltz 8/20/97.
# Name was changed and the "mmt" entry was removed.
class MMTO(Observatory):
	"""
	MMT Observatory
	"""
	def __init__(self):
		self.name      = "MMT Observatory"
		self.longitude = 110.885 * u.degree # West
		self.latitude  = 31.6883333333 * u.degree # North
		self.altitude  = 2600 * u.meter
		self.timezone  = 7 * u.hourangle

class CFHT(Observatory):
	"""
	Canada-France-Hawaii Telescope
	"""
	def __init__(self):
		self.name      = "Canada-France-Hawaii Telescope"
		self.longitude = 155.471666667 * u.degree # West
		self.latitude  = 19.8266666667 * u.degree # North
		self.altitude  = 4215 * u.meter
		self.timezone  = 10 * u.hourangle

class LAPALMA(Observatory):
	"""
	Roque de los Muchachos, La Palma
	"""
	def __init__(self):
		self.name      = "Roque de los Muchachos, La Palma."
		self.longitude = 17.88 * u.degree # West
		self.latitude  = 28.7583333333 * u.degree # North
		self.altitude  = 2327 * u.meter
		self.timezone  = 0 * u.hourangle

class MSO(Observatory):
	"""
	Mt. Stromlo Observatory
	"""
	def __init__(self):
		self.name      = "Mt. Stromlo Observatory"
		self.longitude = 210.975666667 * u.degree # West
		self.latitude  = -34.67935 * u.degree # North
		self.altitude  = 767 * u.meter
		self.timezone  = -10 * u.hourangle

class SSO(Observatory):
	"""
	Siding Spring Observatory
	"""
	def __init__(self):
		self.name      = "Siding Spring Observatory"
		self.longitude = 210.938805556 * u.degree # West
		self.latitude  = -30.7266388889 * u.degree # North
		self.altitude  = 1149 * u.meter
		self.timezone  = -10 * u.hourangle

class AAO(Observatory):
	"""
	Anglo-Australian Observatory
	"""
	def __init__(self):
		self.name      = "Anglo-Australian Observatory"
		self.longitude = 210.933913889 * u.degree # West
		self.latitude  = -30.7229611111 * u.degree # North
		self.altitude  = 1164 * u.meter
		self.timezone  = -10 * u.hourangle

class MCDONALD(Observatory):
	"""
	McDonald Observatory
	"""
	def __init__(self):
		self.name      = "McDonald Observatory"
		self.longitude = 104.0216667 * u.degree # West
		self.latitude  = 30.6716667 * u.degree # North
		self.altitude  = 2075 * u.meter
		self.timezone  = 6 * u.hourangle

class LCO(Observatory):
	"""
	Las Campanas Observatory
	"""
	def __init__(self):
		self.name      = "Las Campanas Observatory"
		self.longitude = 70.7016666667 * u.degree # West
		self.latitude  = -28.9966666667 * u.degree # North
		self.altitude  = 2282 * u.meter
		self.timezone  = 4 * u.hourangle

# Submitted by Alan Koski 1/13/93
class MTBIGELOW(Observatory):
	"""
	Catalina Observatory: 61 inch telescope
	"""
	def __init__(self):
		self.name      = "Catalina Observatory: 61 inch telescope"
		self.longitude = 110.731666667 * u.degree # West
		self.latitude  = 32.4166666667 * u.degree # North
		self.altitude  = 2510. * u.meter
		self.timezone  = 7 * u.hourangle

# Revised by Daniel Durand 2/23/93
class DAO(Observatory):
	"""
	Dominion Astrophysical Observatory
	"""
	def __init__(self):
		self.name      = "Dominion Astrophysical Observatory"
		self.longitude = 123.416666667 * u.degree # West
		self.latitude  = 48.5216666667 * u.degree # North
		self.altitude  = 229 * u.meter
		self.timezone  = 8 * u.hourangle

# Submitted by Patrick Vielle 5/4/93
class SPM(Observatory):
	"""
	Observatorio Astronomico Nacional, San Pedro Martir
	"""
	def __init__(self):
		self.name      = "Observatorio Astronomico Nacional, San Pedro Martir."
		self.longitude = 115.486944444 * u.degree # West
		self.latitude  = 31.0291666667 * u.degree # North
		self.altitude  = 2830 * u.meter
		self.timezone  = 7 * u.hourangle

# Submitted by Patrick Vielle 5/4/93
class TONA(Observatory):
	"""
	Observatorio Astronomico Nacional, Tonantzintla
	"""
	def __init__(self):
		self.name      = "Observatorio Astronomico Nacional, Tonantzintla."
		self.longitude = 98.3138888889 * u.degree # West
		self.latitude  = 19.0327777778 * u.degree # North
		self.timezone  = 8 * u.hourangle

# Submitted by Don Hamilton 8/18/93
class PALOMAR(Observatory):
	"""
	The Hale Telescope
	"""
	def __init__(self):
		self.name      = "The Hale Telescope"
		self.longitude = 116.863 * u.degree # West
		self.latitude  = 33.356 * u.degree # North
		self.altitude  = 1706. * u.meter
		self.timezone  = 8 * u.hourangle

# Submitted by Pat Seitzer 10/31/93
class MDM(Observatory):
	"""
	Michigan-Dartmouth-MIT Observatory
	"""
	def __init__(self):
		self.name      = "Michigan-Dartmouth-MIT Observatory"
		self.longitude = 111.616666667 * u.degree # West
		self.latitude  = 31.95 * u.degree # North
		self.altitude  = 1938.5 * u.meter
		self.timezone  = 7 * u.hourangle

# Submitted by Ignacio Ferrin 9/1/94
class NOV(Observatory):
	"""
	National Observatory of Venezuela
	"""
	def __init__(self):
		self.name      = "National Observatory of Venezuela"
		self.longitude = 70.8666666667 * u.degree # West
		self.latitude  = 8.79 * u.degree # North
		self.altitude  = 3610 * u.meter
		self.timezone  = 4 * u.hourangle

# Submitted by Alan Welty 10/28/94
class BMO(Observatory):
	"""
	Black Moshannon Observatory
	"""
	def __init__(self):
		self.name      = "Black Moshannon Observatory"
		self.longitude = 78.005 * u.degree # West
		self.latitude  = 40.9216666667 * u.degree # North
		self.altitude  = 738 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Biwei JIANG 11/28/95
class BAO(Observatory):
	"""
	Beijing XingLong Observatory
	"""
	def __init__(self):
		self.name      = "Beijing XingLong Observatory"
		self.longitude = 242.425 * u.degree # West
		self.latitude  = 40.3933333333 * u.degree # North
		self.altitude  = 950. * u.meter
		self.timezone  = -8 * u.hourangle

# From Astronomical Almanac 1996
class KECK(Observatory):
	"""
	W. M. Keck Observatory
	"""
	def __init__(self):
		self.name      = "W. M. Keck Observatory"
		self.longitude = 155.478333333 * u.degree # West
		self.latitude  = 19.8283333333 * u.degree # North
		self.altitude  = 4160 * u.meter
		self.timezone  = 10 * u.hourangle

# Submitted by Lina Tomasella 6/11/96:  Padova Astronomical Obs., Asiago, Italy.
class EKAR(Observatory):
	"""
	Mt. Ekar 182 cm. Telescope
	"""
	def __init__(self):
		self.name      = "Mt. Ekar 182 cm. Telescope"
		self.longitude = 348.418866667 * u.degree # West
		self.latitude  = 45.8485888889 * u.degree # North
		self.altitude  = 1413.69 * u.meter
		self.timezone  = -1 * u.hourangle

# Submitted by Michael Ledlow 8/8/96
class APO(Observatory):
	"""
	Apache Point Observatory
	"""
	def __init__(self):
		self.name      = "Apache Point Observatory"
		self.longitude = 105.82 * u.degree # West
		self.latitude  = 32.78 * u.degree # North
		self.altitude  = 2798. * u.meter
		self.timezone  = 7 * u.hourangle

# Submitted by Michael Ledlow 8/8/96
class LOWELL(Observatory):
	"""
	Lowell Observatory
	"""
	def __init__(self):
		self.name      = "Lowell Observatory"
		self.longitude = 111.535 * u.degree # West
		self.latitude  = 35.0966666667 * u.degree # North
		self.altitude  = 2198. * u.meter
		self.timezone  = 7 * u.hourangle

# Submitted by S.G. Bhargavi 8/12/96
class VBO(Observatory):
	"""
	Vainu Bappu Observatory
	"""
	def __init__(self):
		self.name      = "Vainu Bappu Observatory"
		self.longitude = 281.1734 * u.degree # West
		self.latitude  = 12.57666 * u.degree # North
		self.altitude  = 725. * u.meter
		self.timezone  = -5.5 * u.hourangle

# Submitted by S. Giridhar 6/28/03
class IAO(Observatory):
	"""
	Indian Astronomical Observatory, Hanle
	"""
	def __init__(self):
		self.name      = "Indian Astronomical Observatory, Hanle"
		self.longitude = 281.03583 * u.degree # West
		self.latitude  = 32.7794 * u.degree # North
		self.altitude  = 4500 * u.meter
		self.timezone  = -5.5 * u.hourangle

# Submitted by Doug Mink 1/6/97
class FLWO(Observatory):
	"""
	Whipple Observatory
	"""
	def __init__(self):
		self.name      = "Whipple Observatory"
		self.longitude = 110.8775 * u.degree # West
		self.latitude  = 31.6809444444 * u.degree # North
		self.altitude  = 2320 * u.meter
		self.timezone  = 7 * u.hourangle

class FLWO1(Observatory):
	"""
	Whipple Observatory
	"""
	def __init__(self):
		self.name      = "Whipple Observatory"
		self.longitude = 110.8775 * u.degree # West
		self.latitude  = 31.6809444444 * u.degree # North
		self.altitude  = 2320 * u.meter
		self.timezone  = 7 * u.hourangle

# Submitted by Doug Mink 1/6/97
class ORO(Observatory):
	"""
	Oak Ridge Observatory
	"""
	def __init__(self):
		self.name      = "Oak Ridge Observatory"
		self.longitude = 71.5581444444 * u.degree # West
		self.latitude  = 42.5052611111 * u.degree # North
		self.altitude  = 184 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Claudia Vilega Rodriques 12/12/97
class LNA(Observatory):
	"""
	Laboratorio Nacional de Astrofisica - Brazil
	"""
	def __init__(self):
		self.name      = "Laboratorio Nacional de Astrofisica - Brazil"
		self.longitude = 45.5825 * u.degree # West
		self.latitude  = -21.4655555556 * u.degree # North
		self.altitude  = 1864 * u.meter
		self.timezone  = 3 * u.hourangle

# Submitted by John Memzies 12/31/99
class SAAO(Observatory):
	"""
	South African Astronomical Observatory
	"""
	def __init__(self):
		self.name      = "South African Astronomical Observatory"
		self.longitude = 339.189305556 * u.degree # West
		self.latitude  = -31.6205555556 * u.degree # North
		self.altitude  = 1798 * u.meter
		self.timezone  = -2 * u.hourangle

# Submitted by Jorge Federico Gonzalez 12/10/98
class CASLEO(Observatory):
	"""
	Complejo Astronomico El Leoncito, San Juan
	"""
	def __init__(self):
		self.name      = "Complejo Astronomico El Leoncito, San Juan."
		self.longitude = 69.3 * u.degree # West
		self.latitude  = -30.2008333333 * u.degree # North
		self.altitude  = 2552 * u.meter
		self.timezone  = 3 * u.hourangle

# Submitted by Jorge Federico Gonzalez 12/10/98
class BOSQUE(Observatory):
	"""
	Estacion Astrofisica Bosque Alegre, Cordoba
	"""
	def __init__(self):
		self.name      = "Estacion Astrofisica Bosque Alegre, Cordoba."
		self.longitude = 64.5458333333 * u.degree # West
		self.latitude  = -30.4016666667 * u.degree # North
		self.altitude  = 1250 * u.meter
		self.timezone  = 3 * u.hourangle

# Submitted by Ilian Iliev 1/19/99
class ROZHEN(Observatory):
	"""
	National Astronomical Observatory Rozhen - Bulgaria
	"""
	def __init__(self):
		self.name      = "National Astronomical Observatory Rozhen - Bulgaria."
		self.longitude = 335.256111111 * u.degree # West
		self.latitude  = 41.6930555556 * u.degree # North
		self.altitude  = 1759 * u.meter
		self.timezone  = -2 * u.hourangle

# Submitted by Bill Vacca 7/14/99
class IRTF(Observatory):
	"""
	NASA Infrared Telescope Facility
	"""
	def __init__(self):
		self.name      = "NASA Infrared Telescope Facility"
		self.longitude   = 155.471999 * u.degree # West
		self.latitude  = 19.826218 * u.degree # North
		self.altitude  = 4168 * u.meter
		self.timezone  = 10 * u.hourangle

# Submitted by Andy Layden 7/16/99
class BGSUO(Observatory):
	"""
	Bowling Green State Univ Observatory
	"""
	def __init__(self):
		self.name      = "Bowling Green State Univ Observatory."
		self.longitude = 83.6591666667 * u.degree # West
		self.latitude  = 41.3783333333 * u.degree # North
		self.altitude  = 225. * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Oliver-Mark Cordes 8/5/99
class DSAZ(Observatory):
	"""
	Deutsch-Spanisches Observatorium Calar Alto - Spain
	"""
	def __init__(self):
		self.name      = "Deutsch-Spanisches Observatorium Calar Alto - Spain."
		self.longitude = 2.54625 * u.degree # West
		self.latitude  = 37.2236111111 * u.degree # North
		self.altitude  = 2168 * u.meter
		self.timezone  = -1 * u.hourangle

# Submitted by Matilde Fernandez 2/2/99
class CA(Observatory):
	"""
	Calar Alto Observatory
	"""
	def __init__(self):
		self.name      = "Calar Alto Observatory"
		self.longitude = 2.54625 * u.degree # West
		self.latitude  = 37.2236111111 * u.degree # North
		self.altitude  = 2168 * u.meter
		self.timezone  = -1 * u.hourangle

# Submitted by Oliver-Mark Cordes 8/5/99
class HOLI(Observatory):
	"""
	Observatorium Hoher List (Universitaet Bonn) - Germany
	"""
	def __init__(self):
		self.name      = "Observatorium Hoher List (Universitaet Bonn) - Germany."
		self.longitude = 6.85 * u.degree # West
		self.latitude  = 50.16276 * u.degree # North
		self.altitude  = 541 * u.meter
		self.timezone  = -1 * u.hourangle

# Submitted by Steven Majewski 8/27/99
class LMO(Observatory):
	"""
	Leander McCormick Observatory
	"""
	def __init__(self):
		self.name      = "Leander McCormick Observatory"
		self.longitude = 78.5233333333 * u.degree # West
		self.latitude  = 38.0333333333 * u.degree # North
		self.altitude  = 264 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Steven Majewski 8/27/99
class FMO(Observatory):
	"""
	Fan Mountain Observatory
	"""
	def __init__(self):
		self.name      = "Fan Mountain Observatory"
		self.longitude = 78.6933333333 * u.degree # West
		self.latitude  = 37.8783333333 * u.degree # North
		self.altitude  = 566 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Kim K. McLeod 10/13/1999
class WHITIN(Observatory):
	"""
	Whitin Observatory, Wellesley College
	"""
	def __init__(self):
		self.name      = "Whitin Observatory,Wellesley College"
		self.longitude = 71.305833 * u.degree # West
		self.latitude  = 42.295 * u.degree # North
		self.altitude  = 32 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Nuno Peixinho 6/7/2000
# Parameters for the Sierra Nevada Observatory (Spain)
class OSN(Observatory):
	"""
	Observatorio de Sierra Nevada
	"""
	def __init__(self):
		self.name      = "Observatorio de Sierra Nevada"
		self.longitude = 3.38472222222 * u.degree # West
		self.latitude  = 37.0641666667 * u.degree # North
		self.altitude  = 2896 * u.meter
		self.timezone  = -1 * u.hourangle

# Gemini Observatory - Submitted by Inger Jorgensen
#
# Gemini-North
# These values are from the aerial survey made Sept 25, 1996
#    http://www.ifa.hawaii.edu/mko/coordinates.html

class GEMININORTH(Observatory):
	"""
	Gemini North Observatory
	"""
	def __init__(self):
		self.name      = "Gemini North Observatory"
		self.longitude = 155.46904675 * u.degree # West
		self.latitude  = 19.8238015 * u.degree # North
		self.altitude  = 4213.4 * u.meter
		self.timezone  = 10 * u.hourangle

# Gemini-South
# The coordinates are from GPS measurements and the elevation
# is from the Gemini Web pages http://www.gemini.edu/public

# class GEMINI-SOUTH(Observatory):
	# """
	# Gemini South Observatory
	# """
	#def __init__(self):
		# 	self.name      = "Gemini South Observatory"
		# 	self.longitude = 70.7233333333 * u.degree # West
		# 	self.latitude  = -29.7716666667 * u.degree # North
		# 	self.altitude  = 2737. * u.meter
		# 	self.timezone  = 4 * u.hourangle

# Corrected coords from Bryan Miller, 5/18/2006
class GEMINISOUTH(Observatory):
	"""
	Gemini South Observatory
	"""
	def __init__(self):
		self.name      = "Gemini South Observatory"
		self.longitude = 70.7366933333 * u.degree # West
		self.latitude  = -29.75925 * u.degree # North
		self.altitude  = 2722. * u.meter
		self.timezone  = 4 * u.hourangle

# ESO

class LASILLA(Observatory):
	"""
	European Southern Observatory: La Silla
	"""
	def __init__(self):
		self.name      = "European Southern Observatory: La Silla."
		self.longitude = 70.73 * u.degree # West
		self.latitude  = -28.7433333333 * u.degree # North
		self.altitude  = 2347 * u.meter
		self.timezone  = 4 * u.hourangle

class PARANAL(Observatory):
	"""
	European Southern Observatory: Paranal
	"""
	def __init__(self):
		self.name      = "European Southern Observatory: Paranal."
		self.longitude = 70.4033333333 * u.degree # West
		self.latitude  = -23.375 * u.degree # North
		self.altitude  = 2635 * u.meter
		self.timezone  = 4 * u.hourangle

# The following additional entries were suppied by Peter Weilbacher
# weilbach@uni-sw.gwdg.de on 28 Jan 2002 13:12:59 -0700.

class ESONTT(Observatory):
	"""
	European Southern Observatory, NTT, La Silla
	"""
	def __init__(self):
		self.name      = "European Southern Observatory, NTT, La Silla."
		self.longitude = 70.7317422222 * u.degree # West
		self.latitude  = -28.7448777778 * u.degree # North
		self.altitude  = 2375 * u.meter
		self.timezone  = 4 * u.hourangle

class ESO36M(Observatory):
	"""
	European Southern Observatory, 3.6m Telescope, La Silla
	"""
	def __init__(self):
		self.name      = "European Southern Observatory, 3.6m Telescope, La Silla."
		self.longitude = 70.7296127778 * u.degree # West
		self.latitude  = -28.7428294444 * u.degree # North
		self.altitude  = 2400 * u.meter
		self.timezone  = 4 * u.hourangle

class ESOVLT(Observatory):
	"""
	European Southern Observatory, VLT, Paranal
	"""
	def __init__(self):
		self.name      = "European Southern Observatory, VLT, Paranal."
		self.longitude = 70.4022 * u.degree # West
		self.latitude  = -24.6253 * u.degree # North
		self.altitude  = 2648 * u.meter
		self.timezone  = 4 * u.hourangle


# Submited by Giovanni Catanzaro, 7/17/03.
class SLN(Observatory):
	"""
	SLN - Catania Astrophysical Observatory
	"""
	def __init__(self):
		self.name      = "SLN - Catania Astrophysical Observatory."
		self.longitude = 345.026666667 * u.degree # West
		self.latitude  = 37.6916666667 * u.degree # North
		self.altitude  = 1725. * u.meter
		self.timezone  = -1 * u.hourangle


# Submited by Ahmet Devlen, 4/21/04
class EUO(Observatory):
	"""
	Ege University Observatory
	"""
	def __init__(self):
		self.name      = "Ege University Observatory"
		self.longitude = -26.725 * u.degree # West
		self.latitude  = 38.3983333333 * u.degree # North
		self.altitude  = 795. * u.meter
		self.timezone  = 2 * u.hourangle

# Submitted by Zeki Aslan 8/15/05 who said the "tno" entry was wrong.
class TNO(Observatory):
	"""
	Turkiye National Observatory
	"""
	def __init__(self):
		self.name      = "Turkiye National Observatory"
		self.longitude = -29.6469444444 * u.degree # West
		self.latitude  = 36.8244444444 * u.degree # North
		self.altitude  = 2555. * u.meter
		self.timezone  = 2 * u.hourangle

class TUG(Observatory):
	"""
	TUBITAK National Observatory, Turkey
	"""
	def __init__(self):
		self.name      = "TUBITAK National Observatory, Turkey."
		self.longitude = -29.6666666667 * u.degree # West
		self.latitude  = 36.825 * u.degree # North
		self.altitude  = 2547. * u.meter
		self.timezone  = -2 * u.hourangle


# Submited by Ricky Patterson for Vatican Obs. Research Group, 6/15/04
class MGO(Observatory):
	"""
	Mount Graham Observatory
	"""
	def __init__(self):
		self.name      = "Mount Graham Observatory"
		self.longitude = 109.891666667 * u.degree # West
		self.latitude  = 32.7016666667 * u.degree # North
		self.altitude  = 3181 * u.meter
		self.timezone  = 7 * u.hourangle

# Submited by Jeewan C. Bandey 7/28/05
# Changed to W longitude MJF 4/1/06)
# 	(E) longitude = 79.45639
class ARIES(Observatory):
	"""
	Aryabhatta Research Institute of Observational Sciences
	"""
	def __init__(self):
		self.name      = "Aryabhatta Research Institute of Observational Sciences."
		self.longitude = 280.54361 * u.degree # West
		self.latitude  = 29.36 * u.degree # North
		self.altitude  = 1950. * u.meter
		self.timezone  = -5.5 * u.hourangle

# Submitted by Eduardo Fern?ndez Laj?s 10/28/05
class OALP(Observatory):
	"""
	Observatorio Astronomico de La Plata
	"""
	def __init__(self):
		self.name      = "Observatorio Astronomico de La Plata"
		self.longitude = 57.9322995 * u.degree # West
		self.latitude  = -33.0932488889 * u.degree # North
		self.altitude  = 20. * u.meter
		self.timezone  = 3. * u.hourangle

# Submitted by Leslie F. Brown 7/29/06
class OLIN(Observatory):
	"""
	Connecticut College - Olin Observatory
	"""
	def __init__(self):
		self.name      = "Connecticut College - Olin Observatory"
		self.longitude = 72.1052777778 * u.degree # West
		self.latitude  = 41.3788888889 * u.degree # North
		self.altitude  = 85 * u.meter
		self.timezone  = 5 * u.hourangle

# Submitted by Pat van Heerden 11/20/06
class BOYDEN(Observatory):
	"""
	Boyden Observatory
	"""
	def __init__(self):
		self.name      = "Boyden Observatory"
		self.longitude = 332.594444444 * u.degree # West
		self.latitude  = -28.9611111111 * u.degree # North
		self.altitude  = 1387 * u.meter
		self.timezone  = -2 * u.hourangle

# Submitted by Mike Yang 8/19/09
class LULIN(Observatory):
	"""
	Lulin Observatory
	"""
	def __init__(self):
		self.name      = "Lulin Observatory"
		self.longitude = 240.873333333 * u.degree # West
		self.latitude  = 23.4683333333 * u.degree # North
		self.altitude  = 2862. * u.meter
		self.timezone  = -8 * u.hourangle

# Submitted by Mairan Teodoro 1/27/10
class SOAR(Observatory):
	"""
	Southern Astrophysical Research Telescope
	"""
	def __init__(self):
		self.name      = "Southern Astrophysical Research Telescope."
		self.longitude = 70.7337222222 * u.degree # West
		self.latitude  = -29.762 * u.degree # North
		self.altitude  = 2738 * u.meter
		self.timezone  = 4 * u.hourangle

# Submitted from iraf.net 4/12/10
class BAKER(Observatory):
	"""
	Baker Observatory
	"""
	def __init__(self):
		self.name      = "Baker Observatory"
		self.longitude = 93.0417472222 * u.degree # West
		self.latitude  = 37.398625 * u.degree # North
		self.altitude  = 418.2 * u.meter
		self.timezone  = 6 * u.hourangle

# Added MJF 6/1/2010
class HET(Observatory):
	"""
	McDonald Observatory - Hobby-Eberly Telescope
	"""
	def __init__(self):
		self.name      = "McDonald Observatory - Hobby-Eberly Telescope."
		self.longitude = 104.014722222 * u.degree # West
		self.latitude  = 30.68144444 * u.degree # North
		self.altitude  = 2026 * u.meter
		self.timezone  = 6 * u.hourangle

# Submitted by Robert D. Collier 9/1/10
class JCDO(Observatory):
	"""
	Jack C. Davis Observatory, Western Nevada College
	"""
	def __init__(self):
		self.name      = "Jack C. Davis Observatory, Western Nevada College"
		self.longitude = 119.790666667 * u.degree # West
		self.latitude  = 39.1857222222 * u.degree # North
		self.altitude  = 1534 * u.meter
		self.timezone  = 8 * u.hourangle

# Submitted by mas_nomi1711@yahoo.com 3/16/12
class LNO(Observatory):
	"""
	Langkawi National Observatory
	"""
	def __init__(self):
		self.name      = "Langkawi National Observatory"
		self.longitude = 260.218888889 * u.degree # West
		self.latitude  = 6.30694444444 * u.degree # North
		self.altitude  = 111 * u.meter
		self.timezone  = -8 * u.hourangle
