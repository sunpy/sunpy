# -*- coding: utf-8 -*-
# Author:   Michael Malocha
# e-mail:   mmalocha13@gmail.com
# Version:  June 20th, 2013
#

""" This module is to interface the results from a HEK query into a VSO query
    and return the results from the VSO query to the user."""

from __future__ import absolute_import

import json
from sunpy.net import hek

__author__ = "Michael Malocha"

class HEK2VSOTool(object):
    """This tool is to take in a HEK query and return a VSO result"""
    # Still needs the VSO handler

    def __init__(self):
        self.hek_client = hek.HEKClient()

    def query(self, *query):
        """Takes a HEK style query and passes it to a HEKClient instance"""
        results = self.hek_client.query(query)

        with open('h2v-test.log', 'w') as temp_file:
            temp_file.write(json.dumps(results))

        with open('h2v-test.log', 'r') as temp_file:
            results = temp_file.read()

        print json.dumps(json.loads(results), sort_keys=True,
                         indent=4, separators=(',', ': '))

    def testing(self):
        """A docstring"""
        pass

