# -*- coding: utf-8 -*-
# Author:   Michael Malocha
# e-mail:   mmalocha13@gmail.com
# Version:  July 27th, 2013
#

"""
This module is to interface the results from a HEK query into a VSO query
and return the results from the VSO query to the user.
"""

from __future__ import absolute_import

#import json    # Now no longer necessary, though kept till code better tested
import sys
from astropy import units as u
from sunpy.net import hek
from sunpy.net import vso

__author__ = 'Michael Malocha'
__version__ = 'July 17th, 2013'


def wave_unit_catcher(wavelength, units):
    """
    Designed to discover the units of the wavelength passed in and convert
    it into angstroms.

    wavelength: The integer that contains the wavelength value.
    units: The string that contains the units of the wavelength.
    """
    if units == 'nm':
        return u.nm.to(u.angstrom, wavelength)
    elif units == 'mm':
        return u.mm.to(u.angstrom, wavelength)
    elif units == 'cm':
        return u.cm.to(u.angstrom, wavelength)
    elif units == 'm':
        return u.m.to(u.angstrom, wavelength)
    elif units == 'angstrom':
        return u.angstrom.to(u.angstrom, wavelength)
    else:
        print 'Error: Wavelength units not recognized'
        return None


# removed the '*' from results
def translate_results_to_query(results):
    """
    Parses the results from a HEK query and formulates a series of
    VSO queries.

    results: The HEK results from a HEK query to be translated.
    """
    typeList = "<type 'list'>"
    typeHekResponse = "<class 'sunpy.net.hek.hek.Response'>"
    tempType = str(type(results))
    queries = []
    if tempType == typeList:
        for result in results:
            query = vso_attribute_parse(result)
            queries.append(query)
    elif tempType == typeHekResponse:
        query = vso_attribute_parse(results)
        queries.append(query)
    else:
        queries.append(None)
    return queries


def vso_attribute_parse(phrase):
    """
    This is a simple function to parse HEK query result and generate a list
    containing VSO relevant attributes.

    phrase: The single HEK result to be parsed for VSO attribute data.
    """
    query = [vso.attrs.Time(phrase['event_starttime'],
                            phrase['event_endtime']),
             vso.attrs.Source(phrase['obs_observatory']),
             vso.attrs.Instrument(phrase['obs_instrument'])]
    avg_wave_len = wave_unit_catcher(phrase['obs_meanwavel'],
                                     phrase['obs_wavelunit'])
    query.append(vso.attrs.Wave(avg_wave_len, avg_wave_len))
    return query


def progress_bar(phrase, position, total):
    """
    This is a function to print a simple progress bar to the screen

    phrase: The desired phrase to precede the progress bar with each update.
    position: The current location in the data set (will be divided by 'total'
            to determine percent finished).
    total: The total size of the data set
    """
    bar_size = 20
    position = float(format(position, '.1f'))
    fraction = position / total
    percent = str(int(100 * fraction)) + '%'
    place = int(fraction * bar_size)
    pounds = '#' * place
    blanks = ' ' * (20 - place)
    prog_bar = ' [' + pounds + blanks + '] ' + percent
    sys.stdout.write('\r                                                    ')
    sys.stdout.flush()
    sys.stdout.write('\r' + phrase + prog_bar)
    sys.stdout.flush()


class H2VClient(object):
    """
    This tool is to take in a HEK query and return a VSO result
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

    def return_num_records(self):
        """
        Simply returns the number of records from the VSO query
        """
        return self.num_of_records

    def return_vso_results(self):
        """
        Returns the results from the VSO query
        """
        return self.vso_results

    def return_hek_results(self):
        """
        Returns the results from the HEK query
        """
        return self.hek_results

    def full_query(self, client_query, limit=None):
        """
        Takes a list containing a HEK style query, passes it to a HEKClient
         instance, translates it, queries the VSO webservice, then returns
         the VSO results inside a structured list.

         client_query: the list containing the HEK style query.
         limit: an approximate limit to the desired number of VSO results.

         Possible additions:
         - Save results to file
        """
        self.quick_clean()
        sys.stdout.write('\rQuerying HEK webservice...')
        sys.stdout.flush()
        self.hek_results = self.hek_client.query(*client_query)
        return self.translate_and_query(self.hek_results, limit=limit, full_query=True)

    def translate_and_query(self, hek_results, limit=None, full_query=False):
        """
        Takes the results from a HEK query, translates them, then makes a VSO
        query, returning the results in a list organized by their
        corresponding HEK query.

        hek_results: the results from a HEK query in the form of a list.
        limit: an approximate limit to the desired number of VSO results.
        full_query: a simple flag that determines if the method is being
                called from the full_query() method.
        """
        if full_query is False:
            self.quick_clean()
            self.hek_results = hek_results
        vso_query = translate_results_to_query(hek_results)
        result_size = len(vso_query)
        place = 1
        for query in vso_query:
            progress_bar('Querying VSO webservice', place, result_size)
            temp = self.vso_client.query(*query)
            self.vso_results.append(temp)
            self.num_of_records += len(temp)
            if limit is not None:
                if self.num_of_records >= limit:
                    break
            place += 1
        sys.stdout.write('\rDone                                          '
                         '                                                ')
        sys.stdout.flush()
        return self.vso_results

    def quick_clean(self):
        """
        A simple method to quickly sterilize the instance variables.
        """
        self.vso_results = []
        self.num_of_records = 0
