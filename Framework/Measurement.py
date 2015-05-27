# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# See LICENSE (GPLv3)
# slipy/Framework/Measurement.py
"""
Contains the `Measurement` object: A container for a `value` and `error`.

Astropy already has a very useful object, `Quantity` that is expanded into
a `Constant` class. Something with a value, a name, abbreviation, uncertainty,
and a reference. A `Measurement` is nothing more than a `Constant` by
a different name. It functions just like a Quantity/Constant, only we don't
want to be calling it a "constant" and we want to be able to have many of them.
"""

from astropy.units import Quantity

class Measurement(Quantity):
    """
    A container for a `result` and an `error`. This object is meant to be functionally
    equivalent to the astropy.constant.Constant, without the instance checking (we can
    have many `Measurement`s). There are `notes` instead of `references`.
    """
    def __new__(cls, value, error=None, name=None, notes=None):

        instance = super().__new__(cls, value)

        instance.error = error
        instance.name  = name
        instance.notes = notes

        return instance

    # These operations are `broken` because they would otherwise yield a `Measurement`.
    # The right handed operations are not affected, these are all that is needed I think.
    def __truediv__(self, other):
        return Quantity(super().__truediv__(other))

    def __mul__(self, other):
        return Quantity(super().__mul__(other))

    def __add__(self, other):
        return Quantity(super().__add__(other))

    def __sub__(self, other):
        return Quantity(super().__sub__(other))
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __repr__(self):
        return '<' + ' '.join([ label + str(attr) for attr, label in zip(['Measurement',
            self.value * self.unit, self.error, self.name, self.notes],
            ['', '', '| error = ', '| name = ', '| notes = ']) if attr]) + '>'

    def __str__(self):
        attr = [ self.name, self.value*self.unit, self.error, self.notes ]
        name = [' Name  = {}', ' Value = {}', ' Error = {}', ' Notes = {}']
        show = [ a for a in attr if a ]
        return '\n'.join([n for a, n in zip(attr, name) if a]).format(*show)
