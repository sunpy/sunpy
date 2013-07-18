# -*- coding: utf-8 -*-
# Author:   Michael Malocha
# e-mail:   mmalocha13@gmail.com
# Version:  June 20th, 2013
#

"""
This module is to interface the results from a HEK query into a VSO query
and return the results from the VSO query to the user.
"""

from __future__ import absolute_import

import json, sys
from astropy import units as u
from sunpy.net import hek
from sunpy.net import vso

__author__ = 'Michael Malocha'
__version__ = 'July 17th, 2013'


def wave_unit_catcher(wavelength, units):
    """
    Designed to discover the units of the wavelength passed in and convert
    it into angstroms.
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


def translate_results_to_query(*results):
    """
    Parses the results from a HEK query and formulates a series of
    VSO queries.
    """
    results = json.dumps(results)
    results = json.loads(results)
    # The following line is to strip the results from the outer list
    results = results[0]
    queries = []
    for result in results:
        query = [vso.attrs.Time(result['event_starttime'],
                                result['event_endtime']),
                 vso.attrs.Source(result['obs_observatory']),
                 vso.attrs.Instrument(result['obs_instrument'])]
        avg_wave_len = wave_unit_catcher(result['obs_meanwavel'],
                                       result['obs_wavelunit'])
        query.append(vso.attrs.Wave(avg_wave_len, avg_wave_len))
        queries.append(query)
    return queries


def progress_bar(phrase, position, total):
    """
    This is a function to print a simple progress bar to the screen
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

    def query(self, *client_query):
        """
        Takes a HEK style query and passes it to a HEKClient instance.
        Possible additions:
         - Save results to file
        """
        self.num_of_records = 0
        sys.stdout.write('\rQuerying HEK webservice...')
        sys.stdout.flush()
        self.hek_results = self.hek_client.query(*client_query)
        vso_query = translate_results_to_query(self.hek_results)
        result_size = len(vso_query)
        place = 1
        for query in vso_query:
            progress_bar('Querying VSO webservice', place, result_size)
            temp = self.vso_client.query(*query)
            self.vso_results.append(temp)
            self.num_of_records += len(temp)
            place += 1
        sys.stdout.write('\rDone                                          '
                         '                                                ')
        sys.stdout.flush()
        return self.vso_results
