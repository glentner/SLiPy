#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
# slipy/SLiPy/Simbad.py
"""
usage: Simbad.py @Attribute <identifier> [**kwargs]

This module allows the user to query the SIMBAD astronomical database from
inside Python or shell commands/scripts.

The 'Attribute' points to a function within this module and indicates
what is to be run. Execute 'Simbad.py @Attribute help' for usage details of
a specific function. Currently available attributes are: `Position`,
`Distance`, and `Sptype`.

The identifier names can be anything recognized by SIMBAD (e.g., Regulus,
"alpha leo", "HD 121475", "del cyg", etc ...) if the name is two parts make
sure to use quotation to enclose it.

The **kwargs is the conventional reference to Python keyword arguments.
These should be specific to the 'Attribute' being pointed to.
"""
from sys import argv, exit  # , version_info
from urllib.request import urlopen

from astropy import units as u

from .. import SlipyError
from ..Framework.Command import Parse, CommandError
from ..Framework.Options import Options, OptionsError
from ..Framework.Measurement import Measurement


class SimbadError(SlipyError):
    """
    Exception specific to the Simbad module.
    """
pass


def URLEncoded(url):
    """
    URLEncoded( string ):
    Return a string with reserved url characters encoded.
    """
    # check argument type
    if type(url) is not str:
        raise SimbadError('URLEncoded function expects type str!')

    # percent-encoded pair dictionary
    encoded = {
        ' ': '%20',
        '%': '%25',
        '#': '%23',
        '(': '%28',
        ')': '%29',
        '|': '%7c',
        '+': '%2b'
    }

    return ''.join([
        encoded[character] if character in encoded else character
        for character in list(url)])


def Script(identifier, criteria):
    """
    Script( criteria ):

    URL script for the SIMBAD astronomical database with added
    URLEncoded(criteria).
    """
    script = [
        'http://simbad.u-strasbg.fr/simbad/sim-script?',
        'script=format%20object%20%22', URLEncoded(criteria),
        '%22%0a', URLEncoded(identifier)
        ]

    return ''.join(script)


class Query:
    """
    Query( identifier, criteria, **kwargs ):

    Class for querying the SIMBAD astronomical database for 'citeria'
    of 'identifier'.

    kwargs = {
        'parse' : True,  # extract relavent data from SIMBAD return file
        'dtype' : float, # output datatype
    }
    """
    def __init__(self, identifier, criteria, default=float, **kwargs):
        """
        Initiate query to SIMBAD database.
        """
        # check argument types
        if type(identifier) is not str or type(criteria) is not str:
            raise SimbadError('Simbad.Query function expects str'
            'types for arguments.')

        try:
            # keyword argument options for Query
            self.options = Options( kwargs,
                {
                    'parse'  : True    , # parse SIMBAD return file
                    'full'   : False   , # return full line of info
                    'dtype'  : default , # convert return data
                    'is_main': False     # called from Main()
                })

            # assignments
            self.parse   = self.options('parse')
            self.full    = self.options('full')
            self.dtype   = self.options('dtype')
            self.is_main = self.options('is_main')

            # query SIMBAD database
            with urlopen( Script(identifier, criteria) ) as response:
                self.data = str( response.read().decode('utf-8') ).strip()

        except OptionsError as err:
            print('\n --> OptionsError:', err.msg )
            raise SimbadError('Simbad.Query was not constructed '
                'for `{}`'.format(identifier))

        except URLError as error:
            raise SimbadError('Failed to contact SIMBAD database for'
            ' `{}`'.format(identifier) )

        if 'not found' in self.data or 'error' in self.data:
            raise SimbadError('`{}` could not be resolved by SIMBAD.'
                .format(identifier))

        if self.parse:
            # pre-parse operation common to all criteria
            self.data = self.data.split('data')[-1]

    def __call__(self):
        """
        Retrieve data from Query
        """
        return self.data

def Position( identifier, **kwargs ):
    """
    Position( identifier, **kwargs ):

    Handle to the Query class with criteria='%C00(d;C)'.
    """
    query = Query( identifier, '%COO(d;C)', **kwargs )

    if query.full:
        query.data = query.data.split('\n')[-1]

    elif query.parse:
        # extract relavent data
        query.data = query.data.split()[-1].split('+')
        if len( query.data ) == 1:
            # dec had '-' not '+'
            query.data    = query.data[0].split('-')
            query.data[1] = '-' + query.data[1]
        # return formatted data type
        query.data = [ query.dtype(pos) * u.degree for pos in query.data ]

    if query.is_main:
        if query.full or not query.parse:
            print( query() )
        else:
            print('{0:.2f} {1:.2f}'.format(*query()))

    else: return query.data

def Distance( identifier, **kwargs ):
    """
    Distance( identifier, **kwargs ):

    Handle to the Query class with criteria='%PLX'
    """
    query =  Query( identifier, '%PLX', **kwargs )

    if query.full:
        data = query.data.split('\n')[-1]

    elif query.parse:

        # extract relavent data
        data = query.data.split()

        if data[1] == '~':
            # nothing found!
            raise SimbadError('No distance found for `{}`'.format(identifier))
        try:

            # convert milli-arcseconds to parsecs
            result = u.pc / ( query.dtype(data[1]) / 1000.0 )

        except ValueError as err:
            raise SimbadError('Use a numeric type for Simbad.Distance!')

        if data[2][0] == '[' and data[2][-1] == ']':
            # TODO: understand SIMBAD error format
            # an uncertainty was given by SIMBAD
            # uncertainty = u.pc * 0.001 * query.dtype(data[1]) * query.dtype(data[2][1:-1]) / (
            # 	0.001 * query.dtype(data[1]) )**2
            # uncertainty = result * query.dtype(data[2][1:-1])
            uncertainty = None

        else: uncertainty = None

    else: data = query.data

    if query.is_main:
        if query.full or not query.parse:
            print( data )
        else:
            print( '{0:.2f}'.format( data ) )

    elif not query.parse:
        return data

    else: return result

    # Measurement(result, error=uncertainty, name='Distance',
    #     notes='Retrieved from SIMBAD database for `{}`'.format(identifier))

def Sptype(identifier, **kwargs):
    """
    Handle to the Query class with criteria='%SP'.
    """
    query = Query(identifier, '%SP', **kwargs)

    if query.full:
        # return last full line of query
        query.data = query.data.split('\n')[-1]

    elif query.parse:
        # extract relavent data
        query.data = query.data.split()[1]

    if query.is_main:
        print( query() )

    else: return query()

def IDList(identifier, **kwargs):
    """
    Handle to the Query class with criteria='%IDLIST'.
    With `parse` = True, return a list of alternate IDs for
    the `identifier` provided.
    """
    query = Query(identifier, '%IDLIST', **kwargs)

    if query.parse:
        # extract relavent data
        query.data = query.data.split(':')[-1].strip().split('\n')

    if query.is_main:
        for line in query.data:
            print(line)

    else: return query()

def Main( clargs ):
    """
    Main function. See __doc__ for details.
    """

    if len(clargs) < 2:
        # show usage
        print( __doc__ )
        return 0

    # Available functions for execution
    executable = {
            'Distance' : Distance, # search for parsecs
            'Position' : Position, # search for ra, dec
            'Sptype'   : Sptype  , # search for spectral types
            'IDList'   : IDList    # search for IDs
        }


    try:

        # Parse command line arguments
        function, args, kwargs = Parse( clargs[1:] )

        if not args and not kwargs or args[0] == 'help':
            # show function usage
            print( executable[function].__doc__ )
            return 0

        # run execution
        for identifier in args:
            executable[function]( identifier, is_main=True, **kwargs )

        return 0

    except CommandError as err:
        print(' --> CommandError:', err.msg)
        return 1

    except KeyError as key:
        print(' --> {} was not a recognized function.'.format(key))
        return 1

    except SimbadError as err:
        # don't let uncaught self exception pass if from main.
        print(' --> SimbadError:', err.msg)
        return 1

    except Exception as err:
        print(' --> Unrecognized error with query for `{}`'
                .format(args[0]))
        print(' --> Exception: `{}`'.format(err))
        return 1

if __name__ == '__main__':
    # Main return 0 or 1
    exit( Main( argv ) )
