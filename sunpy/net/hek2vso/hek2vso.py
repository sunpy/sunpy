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

import json
from sunpy.net import hek
from sunpy.net import vso

__author__ = "Michael Malocha"


class HEK2VSOTool(object):
    """
    This tool is to take in a HEK query and return a VSO result
    """
    # Still needs the VSO handler

    # Here is some test data
    t_start = '2011/08/09 07:23:56'
    t_end = '2011/08/09 12:40:29'
    t_event = 'FL'

    def __init__(self):
        self.hek_client = hek.HEKClient()
        self.hek_results = ''
        self.vso_client = vso.VSOClient()
        self.vso_results = ''

    def query(self, log=False, log_name=None, *client_query):
        """
        Takes a HEK style query and passes it to a HEKClient instance
        """
        self.hek_results = self.hek_client.query(*client_query)
        vso_query = self.translate_results_to_query(self.hek_results)
        for query in vso_query:
            self.vso_results += self.vso_client.query(query)

        # Saving the results to a local log file. Future changes:
        # - user set path
        # - incorporate timestamp in name
        if log:
            if log_name is not None:
                with open(log_name, 'w') as temp_file:
                    temp_file.write(json.dumps(self.hek_results))
            else:
                with open('h2v-test.log', 'w') as temp_file:
                    temp_file.write(json.dumps(self.hek_results))

        # Temporarily return the initial query results
        return self.vso_results

    def translate_results_to_query(self, results):
        """
        Parses the results from a HEK query and formulates a series of
        VSO queries.
        """
        queries = []
        for result in results:
            query = []
            query.append(vso.attrs.Time(result['event_starttime'],
                                        result['event_endtime']))
            query.append(vso.attrs.Source(result['obs_observatory']))
            query.append(vso.attrs.Instrument(result['obs_instrument']))
            queries.append(query)
        return queries

    def return_hek_results(self, formatted=False, hek_stuff=None):
        """
        A simple method to return the HEK results in either formatted
        or un-formatted JSON form
        """
        if hek_stuff is None:
            results = self.hek_results
        else:
            results = hek_stuff

        if formatted:
            return json.dumps(results, sort_keys=True, indent=4,
                              separators=(',', ': '))
        else:
            return results

    def return_hek_from_file(self, formatted=False):
        """
        A simple testing method to load HEK results from the log file
        """
        temp = formatted
        with open('h2v-test.log', 'r') as temp_file:
            hek_results = json.loads(temp_file.read())
        return self.return_hek_results(formatted=temp, hek_stuff=hek_results)

    def testing(self):
        """
        Testing to see if this class plays nicely with the HEK and
        VSO webservice modules
        """

        # Testing the HEK interface
        client = hek.HEKClient()
        results = client.query(hek.attrs.Time(self.t_start, self.t_end),
                               hek.attrs.EventType(self.t_event))
        if len(results) > 0:
            print "HEK Play-Nice Test: Pass"
        else:
            print "HEK Play-Nice Test: Fail"

        # Testing the VSO interface
        client = vso.VSOClient()
        results = client.query(vso.attrs.Time('2001/1/1', '2001/1/2'),
                               vso.attrs.Instrument('eit'))
        if results.num_records() > 0:
            print "VSO Play-Nice Test: Pass"
        else:
            print "VSO Play-Nice Test: Fail"