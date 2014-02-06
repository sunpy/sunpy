# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 22nd, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

"""
This module is meant to be an interface with the HELIO webservice, providing
it's support through the use of a WSDL file.
"""
from sunpy.net.proxyfix import WellBehavedHttpTransport
from sunpy.net.helio import parser
from sunpy.time import time as T
from suds.client import Client as C
from astropy.io.votable.table import parse_single_table
import io

__author__ = 'Michael Malocha'
__version__ = 'September 22nd, 2013'

# The default wsdl file
DEFAULT_LINK = parser.wsdl_retriever()


def suds_unwrapper(wrapped_data):
    """
    Removes suds wrapping from returned xml data

    When grabbing data via client.last_received() from the suds.client.Client
    module, it returns the xml data in an un-helpful "<s:Envelope>" that needs
    to be removed. This function politely cleans it up.

    Parameters
    ----------
    wrapped_data: str
        Contains the wrapped xml results from a WSDL query

    Returns
    -------
    unwrapped: str
        The xml results with the wrapper removed

    Examples
    --------
    >>> from sunpy.net.helio import hec
    >>> from suds.client import Client
    >>> client = Client(parser.wsdl_retriever())
    >>> client.service.getTableNames()
    >>> temp = client.last_received().str()
    >>> print temp
    <?xml version="1.0" encoding="UTF-8"?>
    <S:Envelope>
       <S:Body>
          <helio:queryResponse>
             <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
                <RESOURCE>
                ...
                </RESOURCE>
             </VOTABLE>
          </helio:queryResponse>
       </S:Body>
    </S:Envelope>
    >>> temp = hec.suds_unwrapper(temp)
    >>> print temp
    <?xml version="1.0" encoding="UTF-8"?>
    <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
        <RESOURCE>
        ...
        </RESOURCE>
     </VOTABLE>
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

    In order to get a VOtable object, it has to be parsed from an xml file or
    file-like object. This function creates a file-like object via the
    StringIO module, writes the xml data to it, then passes the file-like
    object to parse_single_table() from the astropy.io.votable.table module
    and thereby creates a VOtable object.

    Parameters
    ----------
    xml_table: str
        Contains the VOtable style xml data

    Returns
    -------
    votable: astropy.io.votable.tree.Table
        A properly formatted VOtable object

    Examples
    --------
    >>> temp = hec.suds_unwrapper(xml_string)
    >>> type(temp)
    unicode
    >>> temp = hec.votable_handler(temp)
    >>> type(temp)
    astropy.io.votable.tree.Table
    """
    fake_file = io.StringIO()
    fake_file.write(xml_table)
    votable = parse_single_table(fake_file)
    fake_file.close()
    return votable


def time_check_and_convert(time):
    """
    Verifies 'time' as a time object, then makes it VOtable compatible.

    Quickly validates a date-time passed in via a string, then parses out a
    datetime object from the string, which is then converted into the proper
    format to be accepted by the WSDL service methods.
    Returned format: 'yyyy-mm-ddThh:mm:ss'
    Example: '1991-01-03T12:00:00'

    Parameters
    ----------
    time: str
        Contains a singular date-time object within a string.

    Returns
    -------
    converted_time: str or None
        A parsed and converted datetime string.

    Examples
    --------
    >>> time = '2013/01/03'
    >>> hec.time_check_convert(time)
    '2013-01-03T00:00:00'
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
    A client class used to interface with and query HELIO webservices.
    """

    def __init__(self, link=DEFAULT_LINK):
        """
        The constructor; establishes the webservice link for the client

        Initializes the client with a weblink

        Parameters
        ----------
        link: str
            Contains URL to valid WSDL endpoint

        Examples
        --------
        >>> hc = hec.Client()
        """
        self.hec_client = C(link, transport=WellBehavedHttpTransport())

    def time_query(self, start_time, end_time, table=None, max_records=None):
        """
        The simple interface to query the wsdl service.

        Used to utilize the service's TimeQuery() method, this is a simple
        interface between the sunpy module library and the web-service's API.

        Parameters
        ----------
        start_time: str
            The datetime where the query window opens

        end_time: str
            The datetime where the query window closes

        table: str
            The table to query from. If the table is unknown, the user will be
            prompted to pick from a list of tables.

        max_records: int
            The maximum number of desired records.

        Returns
        -------
        results: astropy.io.votable.tree.Table
            Table containing the results from the query

        Examples
        --------
        >>> hc = hec.Client()
        >>> start = '2005/01/03'
        >>> end = '2005/12/03'
        >>> temp = hc.time_query(start, end, max_records=10)
        >>> print temp.array
        [ (31463, '2005-01-03T01:37:36', '2005-01-03T01:37:54', '2005-01-03T01:39:00', 717, 982.0, 113.0, 989.0, 84, 22, 9456, 6, 5010320)
         (31464, '2005-01-03T01:51:36', '2005-01-03T01:59:18', '2005-01-03T02:17:24', 717, 989.0, 117.0, 996.0, 1548, 656, 2286912, 12, 5010301)
         (31465, '2005-01-03T03:26:28', '2005-01-03T03:42:50', '2005-01-03T03:46:04', 717, 994.0, 117.0, 1001.0, 1176, 38, 157800, 6, 5010332)
         (31466, '2005-01-03T03:46:04', '2005-01-03T04:07:10', '2005-01-03T04:07:52', 715, -154.0, 124.0, 198.0, 1308, 1328, 2049360, 12, 5010302)
         (31467, '2005-01-03T05:00:24', '2005-01-03T05:00:30', '2005-01-03T05:19:36', 715, -139.0, 107.0, 176.0, 1152, 224, 894816, 6, 5010313)
         (31468, '2005-01-03T06:40:48', '2005-01-03T06:42:46', '2005-01-03T06:50:12', 717, 990.0, 105.0, 996.0, 564, 23, 50782, 6, 5010314)
         (31469, '2005-01-03T08:27:56', '2005-01-03T08:28:26', '2005-01-03T08:29:08', 717, 971.0, 104.0, 977.0, 72, 36, 11197, 6, 5010334)
         (31470, '2005-01-03T09:31:00', '2005-01-03T09:33:34', '2005-01-03T09:34:52', 717, 960.0, 99.0, 965.0, 232, 108, 56486, 6, 5010322)
         (31471, '2005-01-03T09:34:52', '2005-01-03T09:59:46', '2005-01-03T10:06:04', 717, 994.0, 108.0, 1000.0, 1872, 40, 55920, 6, 5010336)
         (31472, '2005-01-03T11:06:48', '2005-01-03T11:07:18', '2005-01-03T11:15:56', 717, 974.0, 116.0, 981.0, 548, 2160, 2240376, 12, 5010304)]
        """
        if table is None:
            table = self.make_table_list()
        start_time = time_check_and_convert(start_time)
        end_time = time_check_and_convert(end_time)
        if start_time is None or end_time is None:
            return None
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

        Returns the names of all the tables that can be queried via the webservice

        Returns
        -------
        tables.array: numpy.ma.core.MaskedArray
            A VOtable table of available tables names

        Examples
        --------
        >>> hc = hec.Client()
        >>> print hc.get_table_names()
        [('hi_cme_list',) ('cactus_stereoa_cme',) ('aad_gle',)
            ...
         ('wind_sw_crossing_time',) ('ulysses_hxr_flare',)
         ('wind_typeii_soho_cme',)]

        """
        self.hec_client.service.getTableNames()
        tables = self.clean_last_received()
        return tables.array

    def make_table_list(self):
        """
        Creates a list of table names and prompts the user for a choice

        This takes the table of table names from get_table_names(), creates a
        list of the names, sorts them, then presents the tables in a
        convenient menu for the user to choose from. It returns a string
        containing the name of the table that the user picked.

        Returns
        -------
        temp: str
            contains the name of the table that the user picked.

        Examples
        --------
        >>> hc.make_table_list()
          1) aad_gle
          2) aastar_list
          3) apstar_list
          4) bas_magnetic_storms
          ...
        108) wind_waves_type_ii_burst
        109) yohkoh_flare_list
        110) yohkoh_hxr_flare
        111) yohkoh_sxt_trace_list

        Please enter a table number between 1 and 111 ('e' to exit): 108
        'wind_waves_type_ii_burst'
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

        This is a simple function just to clean up code. When called, it grabs
        the "last results", runs suds_unwrapper() on the results, followed by
        votable_handler(), then returns the results to the
        clean_last_received()'s calling method. The "last results" are
        actually tied to the client instance that call the webservice.

        Returns
        -------
        results: astropy.io.votable.tree.Table
            A clean and happy VOtable object

        Examples
        --------
        >>> self.hec_client.service.getTableNames()
        >>> tables = self.clean_last_received()
        """
        results = self.hec_client.last_received().str()
        results = suds_unwrapper(results)
        results = votable_handler(results)
        return results

