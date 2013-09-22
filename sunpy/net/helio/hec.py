# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 15th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

"""
This module is meant to be an interface with the HELIO webservice, providing
it's support through the use of a WSDL file.
"""

from sunpy.net.helio import parser
from suds.client import Client as C
from astropy.io.votable.table import parse_single_table
import io

__author__ = 'Michael Malocha'
__version__ = 'September 15th, 2013'

# The default wsdl file
DEFAULT_LINK = parser.wsdl_retriever()


def suds_unwrapper(wrapped_data):
    """
    Removes suds wrapping from returned xml data

    When grabbing data via client.last_received() from the suds.client.Client
    module, it returns the xml data in an un-helpful "<s:Envelope>" that needs
    to be removed. This function politely cleans it up.
    """
    HEADER = '<?xml version="1.0" encoding="UTF-8"?>'
    CATCH_1 = '<helio:queryResponse>'
    CATCH_2 = '</helio:queryResponse>'
    # Now going to find the locations of the CATCHes in the wrapped_data
    pos_1 = wrapped_data.find(CATCH_1)
    size = len(CATCH_1)
    pos_2 = wrapped_data.find(CATCH_2)
    unwrapped = HEADER + wrapped_data[pos_1 + size:pos_2]
    return unwrapped


def votable_handler(xml_table):
    """
    Returns a VOtable object from a VOtable style xml string
    """
    fake_file = io.StringIO()
    fake_file.write(xml_table)
    votable = parse_single_table(fake_file)
    fake_file.close()
    return votable


class Client(object):
    """
    A simple description of the class
    """
    def __init__(self, link=DEFAULT_LINK):
        """
        The initializing function.
        """
        self.hec_client = C(link)

    def time_query(self, start_time, end_time, table):
        """
        The simple interface to query the wsdl service
        """
        pass

    def get_table_names(self):
        """
        Returns a list of the available tables to query
        """
        self.hec_client.service.getTableNames()
        tables = self.hec_client.last_received().str()
        tables = suds_unwrapper(tables)
        tables = votable_handler(tables)
        return tables.array
        #return self.hec_client.service.getTableNames()




"""
STARTTIME = 'yyyy-mm-ddThh:mm:ss'
ENDTIME = 'yyyy-mm-ddThh:mm:ss'
FROM = table_name
        #Found suffixing 'ivo://helio-vo.eu/hec/table/' in the xml file
        # from the Registry_link

# use sunpy time object for users to enter date for the HEC query
    parse_time - has object to create VO time format

# Revisit and update original co-ordinates converter and re-open pull request
"""
