# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/Data/Atomic.py

"""
Access methods for the atomic data published by Donald C. Morton (2003).
See .Archives.AtomicData.MortonTable
"""

import numpy as np
from astropy import units as u

from .. import SlipyError
from ..Framework.Options import Options, OptionsError
from .Archives import AtomicData

class AtomicError(SlipyError):
    """
    Exception specific for the IonManager class.
    """
    pass

class IonManager:
    """
    Managing class for the atomic data (Morton 2003). See SLiPy.Data.Archives.AtomicData.
    """
    def __init__(self):
        """
        Imports the atomic data and creates the dictionary.
        """
        self.data = AtomicData.MortonTable

        # march through the table and apply units, compute oscillator strengths
        for a, entry in enumerate(self.data):

            AirWave = None if not entry[0] else entry[0] * u.Angstrom
            VacWave = None if not entry[1] else entry[1] * u.Angstrom # shouldn't happen?
            Ion     = entry[2]

            # conversion from 'cm-1' (Mohr & Taylor 2000)
            ELow    = entry[3] * 1.23984186 * u.eV

            # solve for oscillator strengths
            Logwf   = entry[4]
            os      = None if not Logwf else 10**(Logwf) / VacWave.value

            self.data[a] = [AirWave, VacWave, Ion, ELow, Logwf, os]

        # build dictionary by ion name
        self.ions = { ion : [] for ion in set([entry[2] for entry in self.data]) }
        for entry in self.data:
            self.ions[ entry[2] ].append( entry[:2] + entry[3:] )

    def __call__(self, key, **kwargs):
        """
        Retrieve data from the Archives.AtomicData table. If the `key` is a string
        type it is expected to be the name of an ion (e.g., 'C III'). If the `key` is
        a number it is expected to be a wavelength value (if not with units Angstroms are
        implied). The default is Vacuum wavelength, but Air can be specified with the
        keyword argument `wavelength='Air'`.

        If the key was the name of an ion, all the lines for that ion are returned. If the
        key was a wavelength, the closest line in the table to that wavelength is returned.
        You can request a wavelength range by giving the `key` as a tuple of two wavelengths
        specifying the range.

        The results default to the f-value (a.k.a. the oscillator strength) but can be
        changed with the keyword argument `entry`. Options include, `Air`, `Vacuum`, `Ion`,
        `ELow`, `LOGWF`, and `fvalue`.

        The if either a single pair or a list of pairs: the first element of each pair is
        always a wavelength value (in Air if wavelength='Air' or in Vacuum otherwise), the
        second being the entries requested. The wavelength type is always that used for
        the look-up. That is, Vacuum by default, but if `wavelength='Air'` is given, the
        returns will be in Air wavelengths. Be aware that `None` might be returned if
        there is not Air wavelength for the line. Further, all results will be returned
        as `None` where no data is available.
        """
        try:

            options = Options( kwargs, {

                    'wavelength' : 'vacuum', # alternative is `air`
                    'lookup'     : 'fvalue'  # others: `air`, `vacuum`, `ion`, `elow`, `logwf`
                })

            wavelength = options('wavelength')
            lookup     = options('lookup')

        except OptionsError as err:
            print(' --> OptionsError:', err)
            raise AtomicError('Failed keyword assignment from __call__ to IonManager!')

        if isinstance(key, tuple):

            if len(key) != 2:
                raise AtomicError('tuples expected to be length 2 on __call__ to '
                'IonManager!')

            table = self.Between(*key)

        else:

            table = self.__getitem__(key)

        # map entries to index value
        lookup_options = { 'air': 0, 'vacuum': 1, 'ion': 2, 'elow': 3, 'logwf': 4, 'fvalue' : 5 }

        if lookup not in lookup_options:
            raise AtomicError('`{}` is not an available search option!'.format(lookup))

        if wavelength not in ['air', 'vacuum']:
            raise AtomicError('Only `air` and `vacuum` wavelengths are understood!')

        if isinstance(key, str):

            # alter the index dictionary for column changes
            lookup_options = { 'air': 0, 'vacuum': 1, 'elow': 2, 'logwf': 3, 'fvalue' : 4 }

            if lookup == 'ion':
                # `Ion` column won't be present in the returned `table` !!!
                raise AtomicError('You provided the name of an ion but requested the names '
                'of the ions as the return value!')

        if not isinstance(key, tuple) and not isinstance(key, str):
            # assume this was a single wavelength value, `table` is a single line
            return tuple([ table[ lookup_options[wavelength] ], table[ lookup_options[lookup] ] ])

        return [# return pairs of wavelength and `lookup` for each `line` found
                tuple([ line[ lookup_options[wavelength] ], line[ lookup_options[lookup] ] ])
                for line in table
            ]

    def __getitem__(self, index):
        """
        Access methods.

        If the `index` is a string value giving the name of a particular ion, this method
        returns the entries in the data for that ion.

        If the `index` is a numerical value, this method returns the table entry closest
        in wavelength (vacuum).

        If the `index` is a slice object (e.g., `5850:5950`) it returns the segment of the
        table data on that range. A `step` is not acceptable (e.g., `::0.5`).
        """
        if isinstance(index, slice):

            start, stop, step = index.start, index.stop, index.step

            if step: raise AtomicError('You cannot slice the table with a `step` size!')

            if not start:
                # the first vacuum wavelength in the table
                start = self.data[0][1]

            if not stop:
                # the last vacuum wavelength in the table
                stop = self.data[-1][1]

            return self.Between(start, stop)

        elif isinstance(index, str):

            if index not in self.ions:
                raise AtomicError('`{}` is not a recognized or available ion in the '
                'data!'.format(index))

            return self.ions[index]

        else:

            if hasattr(index, 'unit'):
                index = index.to(u.Angstrom).value

            proximity = (np.array([entry[1].value for entry in self.data]) - index)**2

            return self.data[ proximity.argmin() ]

    def Below(self, wavelength):
        """
        Return all table entries below the given wavelength. If units are not given,
        Angstroms are implied.
        """
        if not hasattr(wavelength, 'unit'):
            wavelength *= u.Angstrom

        if wavelength > self.data[-1][1] or wavelength < self.data[0][1]:
            raise AtomicError('Cannot access wavelengths outside the data available!')

        return [ entry for entry in self.data if entry[1] < wavelength ]

    def Above(self, wavelength):
        """
        Return all table entries below the given wavelength. If units are not given,
        Angstroms are implied.
        """
        if not hasattr(wavelength, 'unit'):
            wavelength *= u.Angstrom

        if wavelength > self.data[-1][1] or wavelength < self.data[0][1]:
            raise AtomicError('Cannot access wavelengths outside the data available!')

        return [ entry for entry in self.data if entry[1] > wavelength ]

    def Between(self, wavelengthA, wavelengthB):
        """
        Return all table entries between the given wavelengths. If units are not given,
        Angstroms are implied.
        """
        if not hasattr(wavelengthA, 'unit'):
            wavelengthA *= u.Angstrom

        if not hasattr(wavelengthB, 'unit'):
            wavelengthB *= u.Angstrom

        if wavelengthA > self.data[-1][1] or wavelengthA < self.data[0][1]:
            raise AtomicError('Cannot access wavelengths outside the data available!')

        if wavelengthB > self.data[-1][1] or wavelengthB < self.data[0][1]:
            raise AtomicError('Cannot access wavelengths outside the data available!')

        if wavelengthA > wavelengthB:
            # that doesn't make sense, switch the order
            wavelengthA, WavelengthB = WavelengthB, WavelengthA

        return [ entry for entry in self.Below(wavelengthB) if entry[1] > wavelengthA ]


# create a static instance of the IonManager class.
IonSearch = IonManager()

class Ion():
    """
    An object for declaring atomic ions.

    An `Ion` should have a name, wavelength (air or vacuum), osciallator strength,
    and a transition probability (Einstein coefficient).

    Much of this data is available from the ..Data.Archives.AtomicData module searchable
    with the ..Data.Atomic.IonManager class. In the future, we hope to impliment an
    internal database of transition probabilities for each ion in the AtomicData module
    from the Morton 2003 publication. Currently, the user must provide this data.

    The NIST Atomic Spectra Database Lines Data is a good source for atomic data values:
        http://physics.nist.gov/PhysRefData/ASD/lines_form.html
    """
    def __init__(self, name=None, wavelength=None, fvalue=None, A=None, **kwargs):
        """
        Create a new `Ion`. If no arguments are given the state remains uninitialized.
        The `name` (e.g., 'Ca III') is used to connect with the data in the
        ..Data.Archives.AtomicData module via the IonManager class. The `wavelength` need
        not necessarily be the exact value of the line; the line closest to that given
        for the `name`d ion is used. This is by default the wavelength in vacuum, but can
        be the wavelength in air if the keyword argument `medium='air'` is given.
        The created member attribute `wavelength` and `fvalue` are automatically assigned
        based on this procedure, but can be explicitely assigned if the `fvalue` is
        provided directly. The transition probability (Einstein coefficient) `A` simply
        attached as `A`. If not units are given for `A`, `s-1` is assigned.
        """
        try:

            options = Options( kwargs, {

                    'medium' : 'vacuum' # alternative is 'air'
                })

            self.medium = options('medium')

            if self.medium not in ['vacuum', 'air']:
                raise AtomicError('From Ion.__init__(), the only allowed values for `medium` '
                'are `vacuum` and `air`!')

        except OptionsError as err:
            print(' --> OptionsError:', err)
            raise AtomicError('From Ion.__init__(), that was not an exceptable keyword assignment '
            'for an `Ion`!')

        if not name or not wavelength or fvalue:
            # simply assign whatever is provided.
            self.name       = name
            self.wavelength = wavelength
            self.fvalue     = fvalue
            self.A          = A

        else:

            if not hasattr(wavelength, 'unit'):
                wavelength *= u.Angstrom

            # lookup the ion using IonManager (already created as IonSearch)
            found_ions  = IonSearch(name, wavelength=self.medium)
            wavelengths = np.array([(wavelength - entry[0]).value for entry in found_ions])**2
            closest_ion = wavelengths.argmin()

            # assign the closest wavelength to that avaiable
            self.wavelength = found_ions[ closest_ion ][0].to(wavelength.unit)
            self.fvalue     = found_ions[ closest_ion ][1] * u.dimensionless_unscaled

            self.name = name

            if A and not hasattr(A, 'unit'):
                A /= u.s

            self.A    = A

    def __str__(self):
        """
        Show the available information for this Ion.
        """
        return ('Ion:        {}\n'
                'Wavelength: {}\n'
                'fvalue:     {}\n'
                'A:          {}\n'.format(self.name, self.wavelength, self.fvalue, self.A))

    def __repr__(self):
        """
        The same as the `str` representation.
        """
        return str(self)
