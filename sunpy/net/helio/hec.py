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
from sunpy.time import time as T
from suds.client import Client as C
from astropy.io.votable.table import parse_single_table
import io
import datetime

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


def time_check_and_convert(time):
    """
    Verifies 'time' as a time object, then makes it VOtable compatible
    """
    if T.is_time(time):
        time_object = T.parse_time(time)
        converted_time = time_object.strftime("%Y-%m-%d") + "T" + \
                         time_object.strftime("%T")
        return converted_time
    else:
        print "Not a valid datetime"
        return None


class Client(object):
    """
    A simple description of the class
    """

    def __init__(self, link=DEFAULT_LINK):
        """
        The initializing function.
        """
        self.hec_client = C(link)

    def time_query(self, start_time, end_time, table=None, max_records=None):
        """
        The simple interface to query the wsdl service
        """
        if table is None:
            table = self.make_table_list()
        if max_records is not None:
            self.hec_client.service.TimeQuery(STARTTIME=start_time,
                                              ENDTIME=end_time,
                                              FROM=table,
                                              MAXRECORDS=max_records)
        else:
            self.hec_client.service.TimeQuery(STARTTIME=start_time,
                                              ENDTIME=end_time,
                                              FROM=table)
        results = self.clean_last_received()
        return results

    def get_table_names(self):
        """
        Returns a list of the available tables to query
        """
        self.hec_client.service.getTableNames()
        tables = self.clean_last_received()
        return tables.array

    def make_table_list(self):
        """
        Creates a list of table names and prompts the user for a choice
        """
        table_list = []
        temp = None
        tables = self.get_table_names()
        for i in tables:
            temp = str(i)[2:-3]
            if len(temp) > 0:
                table_list.append(temp)
        counter = 1
        table_list.sort()
        for i in table_list:
            temp = '  '
            if 9 < counter < 100:
                temp = ' '
            elif 99 < counter:
                temp = ''
            temp += str(counter) + ") " + i
            print temp
            counter += 1
        while True:
            input = raw_input("\nPlease enter a table number between 1 and %i "
                              "('e' to exit): " % len(table_list))
            if input.lower() == "e" or input.lower() == "exit":
                temp = None
                break
            temp = [int(s) for s in input.split() if s.isdigit()]
            temp = temp[0] - 1
            if temp in range(0, len(table_list)):
                temp = table_list[temp]
                break
            else:
                print "Choice outside of bounds"
        return temp

    def clean_last_received(self):
        """
        A method to return clean VOtable objects for the client.
        """
        results = self.hec_client.last_received().str()
        results = suds_unwrapper(results)
        results = votable_handler(results)
        return results


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
