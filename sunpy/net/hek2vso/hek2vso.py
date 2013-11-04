# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  August 10th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#
#pylint: disable=W0142

"""
This module is to interface the results from a HEK query into a VSO query
and return the results from the VSO query to the user.
"""

from __future__ import absolute_import

import sys
from astropy import units
from sunpy.net import hek
from sunpy.net import vso

__author__ = 'Michael Malocha'
__version__ = 'Aug 10th, 2013'


def wave_unit_catcher(wavelength, wave_units):
    """
    Catch and convert wavelength to angstroms.

    Designed to discover the units of the wavelength passed in and convert
    it into angstroms. Returns an integer or None.

    Parameters
    ----------
    wavelength : int
        Wavelength value.
    units : str
        Units of the wavelength.

    Examples
    --------
    >>> wave_unit_catcher(2.11e-06, 'cm')
    210.99999999999997

    >>> wave_unit_catcher(9.4e-07, 'cm')
    93.99999999999999

    >>> wave_unit_catcher(5e-08, 'mm')
    0.4999999999999999
    """
    try:
        converted_value = getattr(units, wave_units).to(units.angstrom,
                                                        wavelength)
    except AttributeError:
        raise AttributeError("'%s' is not a supported unit" % wave_units)
    return converted_value


def translate_results_to_query(results):
    """
    Formulate VSO queries from HEK results.

    Take the results from a HEK query either in the form of a single HEK
    response or a list containing multiple HEK responses then translates
    them into a VSO compatible query.

    Parameters
    ----------
    results : sunpy.net.hek.hek.Response or list of sunpy.net.hek.hek.Response
        The HEK results from a HEK query to be translated.

    Examples
    --------
    >>> h = hek.HEKClient()
    >>> h2v = H2VClient()
    >>> q = h.query(hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL'))
    >>> len(q)
    19

    >>> translate_results_to_query(q[0]) # doctest: +ELLIPSIS
    [[<Time(datetime.datetime(2011, 8, 8, 1, 30, 4), datetime.datetime(2011, 8, 10, 0, 0, 4), None)>, <Source(u'SDO')>, <Instrument(u'AIA')>, <sunpy.net.vso.attrs.Wave at 0x...>]]

    >>> translate_results_to_query(q) # doctest: +ELLIPSIS
    [[<Time(datetime.datetime(2011, 8, 8, 1, 30, 4), datetime.datetime(2011, 8, 10, 0, 0, 4), None)>, <Source(u'SDO')>, <Instrument(u'AIA')>, <sunpy.net.vso.attrs.Wave at 0x...>],
    ...
    [<Time(datetime.datetime(2011, 8, 9, 8, 1, 21), datetime.datetime(2011, 8, 9, 8, 16, 45), None)>, <Source(u'SDO')>, <Instrument(u'AIA')>, <sunpy.net.vso.attrs.Wave at 0x...>]]
    """
    queries = []
    if type(results) is list:
        for result in results:
            query = vso_attribute_parse(result)
            queries.append(query)
    else:
        query = vso_attribute_parse(results)
        queries.append(query)
    return queries


def vso_attribute_parse(phrase):
    """
    Parses VSO attributes from a HEK result.

    This is a simple function to parse HEK query result and generate a list
    containing VSO relevant attributes.

    Parameters
    ----------
    phrase: dictionary containing a sunpy.net.hek.hek.Response
        The single HEK result to be parsed for VSO attribute data.

    Examples
    --------
    >>> h = hek.HEKClient()
    >>> h2v = H2VClient()
    >>> q = h.query(hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL'))
    >>> len(q)
    19

    >>> vso_attribute_parse(q[9])
    [<Time(datetime.datetime(2011, 8, 9, 7, 22, 38), datetime.datetime(2011, 8, 9, 8, 32, 2), None)>,
    <Source(u'SDO')>,
    <Instrument(u'AIA')>,
    <sunpy.net.vso.attrs.Wave at 0x10628f950>]
    """
    try:
        query = [vso.attrs.Time(phrase['event_starttime'],
                                phrase['event_endtime']),
                 vso.attrs.Source(phrase['obs_observatory']),
                 vso.attrs.Instrument(phrase['obs_instrument'])]
        avg_wave_len = wave_unit_catcher(phrase['obs_meanwavel'],
                                         phrase['obs_wavelunit'])
        query.append(vso.attrs.Wave(avg_wave_len, avg_wave_len))
    except KeyError, TypeError:
        raise TypeError("'%s' is an improper data type" % type(phrase))
    return query


def progress_bar(phrase, position, total, bar_size=20):
    """
    Prints a simple progress bar to the screen

    Parameters
    ----------
    phrase: str
        The desired phrase to precede the progress bar with each update.
    position: int
        The current location in the data set (will be divided by 'total'
        to determine percent finished).
    total: int
        The total size of the data set
    bar_size: int
        Size of the bar, defaults to 20 characters in length

    Examples
    --------
    >>> progress_bar("Progress:", 10, 100)
    Progress: [##                  ] 10%

    >>> progress_bar("Progress:", 35, 100)
    Progress: [#######             ] 35%

    >>> progress_bar("Progress:", 35, 83)
    Progress: [########            ] 42%
    """
    position = float(format(position, '.1f'))
    fraction = position / total
    percent = str(int(100 * fraction)) + '%'
    place = int(fraction * bar_size)
    pounds = '#' * place
    blanks = ' ' * (20 - place)
    prog_bar = ' [' + pounds + blanks + '] ' + percent
    sys.stdout.write('\r' + ' ' * 52)
    sys.stdout.flush()
    sys.stdout.write('\r' + phrase + prog_bar)
    sys.stdout.flush()


class H2VClient(object):
    """
    Class to handle HEK to VSO translations

    Though the single step functions exists outside this class where
    translation is also possible, this class provides a framework where
    all the necessary functionality is easily accessed, along with a few
    additional and helpful methods.

    Examples
    --------
    >>> from sunpy.net import hek2vso
    >>> h2v = hek2vso.H2VClient()
    """
    # Here is some test data
    t_start = '2011/08/09 07:23:56'
    t_end = '2011/08/09 12:40:29'
    t_event = 'FL'

    def __init__(self):
        self.hek_client = hek.HEKClient()
        self.hek_results = ''
        self.vso_client = vso.VSOClient()
        self.vso_results = []
        self.num_of_records = 0

    def full_query(self, client_query, limit=None):
        """
        An encompassing method that takes a HEK query and returns a VSO result

        Takes a list containing a HEK style query, passes it to a HEKClient
        instance, translates it, queries the VSO webservice, then returns
        the VSO results inside a structured list.

        Parameters
        ----------
        client_query: list
            The list containing the HEK style query.
        limit: int
            An approximate limit to the desired number of VSO results.

        Examples
        --------
        >>> from sunpy.net import hek, hek2vso
        >>> h2v = H2VClient()
        >>> q = h2v.full_query((hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL')))
        """
        self.quick_clean()
        sys.stdout.write('\rQuerying HEK webservice...')
        sys.stdout.flush()
        self.hek_results = self.hek_client.query(*client_query)
        return self.translate_and_query(self.hek_results,
                                        limit=limit, full_query=True)

    def translate_and_query(self, hek_results, limit=None,
                            full_query=False, use_progress_bar=True):
        """
        Translates HEK results, makes a VSO query, then returns the results.

        Takes the results from a HEK query, translates them, then makes a VSO
        query, returning the results in a list organized by their
        corresponding HEK query.

        Parameters
        ----------
        hek_results: sunpy.net.hek.hek.Response or list of sunpy.-.-.-.Response
            The results from a HEK query in the form of a list.
        limit: int
            An approximate limit to the desired number of VSO results.
        full_query: Boolean
            A simple flag that determines if the method is being
            called from the full_query() method.
        use_progress_bar: Boolean
            A flag to turn off the progress bar, defaults to "on"

        Examples
        --------
        >>> from sunpy.net import hek, hek2vso
        >>> h = hek.HEKClient()
        >>> tstart = '2011/08/09 07:23:56'
        >>> t_end = '2011/08/09 12:40:29'
        >>> t_event = 'FL'
        >>> q = h.query(hek.attrs.Time(tstart, tend), hek.attts.EventType(event_type))
        >>> h2v = hek2vso.H2VClient()
        >>> res = h2v.translate_and_query(q)
        """
        if full_query is False:
            self.quick_clean()
            self.hek_results = hek_results
        vso_query = translate_results_to_query(hek_results)
        result_size = len(vso_query)
        place = 1
        for query in vso_query:
            if use_progress_bar:
                progress_bar('Querying VSO webservice', place, result_size)
            temp = self.vso_client.query(*query)
            self.vso_results.append(temp)
            self.num_of_records += len(temp)
            if limit is not None:
                if self.num_of_records >= limit:
                    break
            place += 1
        if use_progress_bar:
            sys.stdout.write('\rDone                                         '
                             '                                              ')
        sys.stdout.flush()
        return self.vso_results

    def quick_clean(self):
        """
        A simple method to quickly sterilize the instance variables.

        Used to bleach local variables before a new query is made. Not
        intended to be run by user.
        """
        self.vso_results = []
        self.num_of_records = 0
