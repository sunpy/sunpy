"""
Access the Helio Event Catalogue
"""
from sunpy.net.proxyfix import WellBehavedHttpTransport
from sunpy.net.helio import parser
from sunpy.time import parse_time
from suds.client import Client as C
import suds
from astropy.io.votable.table import parse_single_table
import io

__author__ = 'Michael Malocha'
__version__ = 'September 22nd, 2013'

__all__ = ['HECClient']


        

def suds_unwrapper(wrapped_data):
    """
    Removes suds wrapping from returned xml data

    When grabbing data via votable_interceptor.last_payload from the suds.client.Client
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
    >>> from sunpy.net.helio import hec  Todo: Fix this example!
    >>> from suds.client import Client
    >>> from sunpy.net.proxyfix import WellBehavedHttpTransport
    >>> votable_interceptor = hec.VotableInterceptor()
    >>> client = Client(hec.parser.wsdl_retriever(), plugins=[self.votable_interceptor], transport=WellBehavedHttpTransport())
    >>> client.service.getTableNames()
    >>> temp = client.last_received().str()
    >>> print temp
    <?xml version="1.0" encoding="UTF-8"?>
    <S:Envelope ..... >
       <S:Body>
          <helio:queryResponse ... >
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
    HEADER = '<?xml version="1.0" encoding="UTF-8"?>\n'
    CATCH_1 = '<VOTABLE'
    CATCH_2 = '</VOTABLE>\n'
    # Now going to find the locations of the CATCHes in the wrapped_data
    pos_1 = wrapped_data.find(CATCH_1)
    pos_2 = wrapped_data.find(CATCH_2)
    unwrapped = HEADER + wrapped_data[pos_1:pos_2] + CATCH_2
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


class VotableInterceptor(suds.plugin.MessagePlugin):
    '''
    Adapted example from http://stackoverflow.com/questions/15259929/configure-suds-to-use-custom-response-xml-parser-for-big-response-payloads
    '''
    def __init__(self, *args, **kwargs):
        self.last_payload = None

    def received(self, context):
        #recieved xml as a string
        self.last_payload = unicode(suds_unwrapper(context.reply))
        #clean up reply to prevent parsing
        context.reply = ""
        return context


class HECClient(object):
    """
    A client class used to interface with and query HELIO webservices.
    """

    def __init__(self, link=None):
        """
        The constructor; establishes the webservice link for the client

        Initializes the client with a weblink

        Parameters
        ----------
        link: str
            Contains URL to valid WSDL endpoint

        Examples
        --------
        >>> hc = hec.HECClient()
        """
        if link is None:
            # The default wsdl file
            link = parser.wsdl_retriever()

        self.votable_interceptor = VotableInterceptor()
        self.hec_client = C(link, plugins=[self.votable_interceptor], transport=WellBehavedHttpTransport())

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
        >>> hc = hec.HECClient()
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
        while table is None:
            table = self.make_table_list()
        start_time = parse_time(start_time)
        end_time = parse_time(end_time)
        self.hec_client.service.TimeQuery(STARTTIME=start_time.isoformat(),
                                          ENDTIME=end_time.isoformat(),
                                          FROM=table,
                                          MAXRECORDS=max_records)
        results = votable_handler(self.votable_interceptor.last_payload)
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
        >>> hc = hec.HECClient()
        >>> print hc.get_table_names()
        [('hi_cme_list',) ('cactus_stereoa_cme',) ('aad_gle',)
            ...
         ('wind_sw_crossing_time',) ('ulysses_hxr_flare',)
         ('wind_typeii_soho_cme',)]

        """
        self.hec_client.service.getTableNames()
        tables = votable_handler(self.votable_interceptor.last_payload)
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
        tables = self.get_table_names()
        for i in tables:
            table = i[0]
            if len(table) > 0:
                table_list.append(table)
        table_list.sort()
        for index, table in enumerate(table_list):
            print ('{number:3d}) {table}'.format(number = index + 1, table = table))
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

