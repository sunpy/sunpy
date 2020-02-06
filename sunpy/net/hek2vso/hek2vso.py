# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  August 10th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#
#pylint: disable=W0142

"""
This module translates the results of a HEK query into a VSO query
and returns the results from the VSO query to the user.
"""

import sys
from astropy import units
from astropy.table import Table

from tqdm import tqdm

from sunpy.net import hek
from sunpy.net import vso
from sunpy.net import attrs as a

__author__ = 'Michael Malocha'
__version__ = 'Aug 10th, 2013'

__all__ = ['translate_results_to_query', 'vso_attribute_parse', 'H2VClient']


def translate_results_to_query(results):
    """
    Formulate VSO queries from HEK results.

    Take the results from a HEK query either in the form of a single HEK
    response or a list containing multiple HEK responses then translates
    them into a VSO compatible query.

    Parameters
    ----------
    results : `sunpy.net.hek.hek.HEKRow` or `sunpy.net.hek.hek.HEKTable`
        The HEK results from a HEK query to be translated.

    Examples
    --------
    >>> from sunpy.net.hek import hek, HEKClient
    >>> from sunpy.net.hek2vso import hek2vso, H2VClient
    >>> h = HEKClient()  # doctest: +REMOTE_DATA
    >>> h2v = H2VClient()  # doctest: +REMOTE_DATA
    >>> q = h.search(hek.attrs.Time('2011/08/09 07:23:56',
    ...             '2011/08/09 12:40:29'), hek.attrs.EventType('FL'))  # doctest: +REMOTE_DATA
    >>> len(q)  # doctest: +REMOTE_DATA
    19

    >>> hek2vso.translate_results_to_query(q[0])  # doctest: +REMOTE_DATA
    [[<Time(<Time object: scale='utc' format='isot' value=2011-08-08T01:30:04.000>, <Time object: scale='utc' format='isot' value=2011-08-10T00:00:04.000>, None)>, <Source('SDO')>, <Instrument('HMI')>, <Wavelength(6172.999999999998, 6172.999999999998, 'Angstrom')>]]
    >>> hek2vso.translate_results_to_query(q)   # doctest: +SKIP
    [[<Time(<Time object: scale='utc' format='isot' value=2011-08-08T01:30:04.000>, <Time object: scale='utc' format='isot' value=2011-08-10T00:00:04.000>, None)>, <Source('SDO')>, <Instrument('HMI')>, <Wavelength(6172.999999999998, 6172.999999999998, 'Angstrom')>], ..., [<Time(<Time object: scale='utc' format='isot' value=2011-08-09T08:01:21.000>, <Time object: scale='utc' format='isot' value=2011-08-09T08:16:45.000>, None)>, <Source('SDO')>, <Instrument('AIA')>, <Wavelength(303.99999999999994, 303.99999999999994, 'Angstrom')>]]
    """
    queries = []
    if isinstance(results, Table):
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
    phrase : `sunpy.net.hek.hek.HEKRow`.
        The single HEK result to be parsed for VSO attribute data.

    Examples
    --------
    >>> from sunpy.net.hek import hek, HEKClient
    >>> from sunpy.net.hek2vso import hek2vso, H2VClient
    >>> h = HEKClient()  # doctest: +REMOTE_DATA
    >>> h2v = H2VClient()  # doctest: +REMOTE_DATA
    >>> q = h.search(hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL'))  # doctest: +REMOTE_DATA
    >>> len(q)  # doctest: +REMOTE_DATA
    19

    >>> hek2vso.vso_attribute_parse(q[9])  # doctest: +REMOTE_DATA
    [<Time(<Time object: scale='utc' format='isot' value=2011-08-09T07:22:38.000>, <Time object: scale='utc' format='isot' value=2011-08-09T08:32:02.000>, None)>, <Source('SDO')>, <Instrument('AIA')>, <Wavelength(210.99999999999997, 210.99999999999997, 'Angstrom')>]
    """
    try:
        query = [a.Time(phrase['event_starttime'],
                                phrase['event_endtime']),
                 a.vso.Source(phrase['obs_observatory']),
                 a.Instrument(phrase['obs_instrument'])]
        avg_wave_len = phrase['obs_meanwavel'] * units.Unit(phrase['obs_wavelunit'])
        query.append(a.Wavelength(avg_wave_len, avg_wave_len))
    except (KeyError, TypeError):
        raise TypeError("'{dtype!s}' is an improper data type".format(dtype=type(phrase)))
    return query


class H2VClient:
    """
    Class to handle HEK to VSO translations

    Though the single step functions exists outside this class where
    translation is also possible, this class provides a framework where
    all the necessary functionality is easily accessed, along with a few
    additional and helpful methods.

    Examples
    --------
    >>> from sunpy.net.hek import hek
    >>> from sunpy.net import hek2vso
    >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA
    """

    def __init__(self):
        self.hek_client = hek.HEKClient()
        self.hek_results = ''
        self.vso_client = vso.VSOClient()
        self.vso_results = []
        self.num_of_records = 0

    def full_query(self, client_query, limit=None, progress=False):
        """
        An encompassing method that takes a HEK query and returns a VSO result

        Takes a list containing a HEK style query, passes it to a HEKClient
        instance, translates it, queries the VSO webservice, then returns
        the VSO results inside a structured list.

        Parameters
        ----------
        client_query : `list`
            The list containing the HEK style query.
        limit : `int`
            An approximate limit to the desired number of VSO results.

        Examples
        --------
        >>> from sunpy.net import hek, hek2vso
        >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA
        >>> q = h2v.full_query((hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'),
        ...                    hek.attrs.EventType('FL')))  # doctest: +REMOTE_DATA
        """
        self._quick_clean()
        if progress:
            sys.stdout.write('\rQuerying HEK webservice...')
            sys.stdout.flush()
        self.hek_results = self.hek_client.search(*client_query)
        self._quick_clean()
        return self.translate_and_query(self.hek_results,
                                        limit=limit, progress=progress)

    def translate_and_query(self, hek_results, limit=None, progress=False):
        """
        Translates HEK results, makes a VSO query, then returns the results.

        Takes the results from a HEK query, translates them, then makes a VSO
        query, returning the results in a list organized by their
        corresponding HEK query.

        Parameters
        ----------
        hek_results : `sunpy.net.hek.hek.HEKRow` or `sunpy.net.hek.hek.HEKTable`
            The results from a HEK query in the form of a list.
        limit : int
            An approximate limit to the desired number of VSO results.
        progress : Boolean
            A flag to turn off the progress bar, defaults to "off"

        Examples
        --------
        >>> from sunpy.net import hek, hek2vso
        >>> h = hek.HEKClient()  # doctest: +REMOTE_DATA
        >>> tstart = '2011/08/09 07:23:56'
        >>> tend = '2011/08/09 12:40:29'
        >>> event_type = 'FL'
        >>> q = h.search(hek.attrs.Time(tstart, tend), hek.attrs.EventType(event_type))  # doctest: +REMOTE_DATA
        >>> h2v = hek2vso.H2VClient()  # doctest: +REMOTE_DATA
        >>> res = h2v.translate_and_query(q)  # doctest: +REMOTE_DATA
        """
        vso_query = translate_results_to_query(hek_results)

        for query in tqdm(vso_query, unit="records"):
            temp = self.vso_client.search(*query)
            self.vso_results.append(temp)
            self.num_of_records += len(temp)
            if limit is not None:
                if self.num_of_records >= limit:
                    break

        return self.vso_results

    def _quick_clean(self):
        """
        A simple method to quickly sterilize the instance variables.

        Used to bleach local variables before a new query is made. Not
        intended to be run by user.
        """
        self.vso_results = []
        self.num_of_records = 0
