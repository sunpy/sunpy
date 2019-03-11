"""
Access the Helio Event Catalogue
"""
import io

import zeep
from lxml import etree

from astropy.io.votable.table import parse_single_table

from sunpy.time import parse_time
from sunpy.net.helio import parser

__all__ = ['HECClient']


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
    xml_table : `bytes`
        Contains the VOtable style xml data

    Returns
    -------
    votable : `astropy.io.votable.tree.Table`
        A properly formatted VOtable object

    """
    fake_file = io.BytesIO()
    fake_file.write(xml_table)
    votable = parse_single_table(fake_file)
    fake_file.close()
    return votable


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
        link : str
            Contains URL to valid WSDL endpoint

        Examples
        --------
        >>> from sunpy.net.helio import hec
        >>> hc = hec.HECClient()  # doctest: +REMOTE_DATA

        """
        if link is None:
            # The default wsdl file
            link = parser.wsdl_retriever()

        self.hec_client = zeep.Client(link)

    def time_query(self, start_time, end_time, table=None, max_records=None):
        """
        The simple interface to query the wsdl service.

        Used to utilize the service's TimeQuery() method, this is a simple
        interface between the sunpy module library and the web-service's API.

        Parameters
        ----------
        start_time : str, `~sunpy.time.parse_time` parsable objects
            The time where the query window opens

        end_time : str, `~sunpy.time.parse_time` parsable objects
            The time where the query window closes

        table : bytes
            The table to query from. If the table is unknown, the user will be
            prompted to pick from a list of tables.

        max_records: int
            The maximum number of desired records.

        Returns
        -------
        results: `astropy.io.votable.tree.Table`
            Table containing the results from the query

        Examples
        --------
        >>> from sunpy.net.helio import hec
        >>> hc = hec.HECClient()  # doctest: +REMOTE_DATA
        >>> start = '2005/01/03'
        >>> end = '2005/12/03'
        >>> temp = hc.time_query(start, end, max_records=10)   # doctest: +REMOTE_DATA +SKIP

        """
        while table is None:
            table = self.select_table()
        start_time = parse_time(start_time)
        end_time = parse_time(end_time)
        results = self.hec_client.service.TimeQuery(STARTTIME=start_time.isot,
                                                    ENDTIME=end_time.isot,
                                                    FROM=table,
                                                    MAXRECORDS=max_records)
        results = votable_handler(etree.tostring(results))
        return results

    def get_table_names(self):
        """
        Returns a list of the available tables to query.

        Returns the names of all the tables that can be queried via the
        webservice.

        Returns
        -------
        tables.array: `numpy.ma.core.MaskedArray`
            A VOtable table of available tables names

        Examples
        --------
        >>> from sunpy.net.helio import hec
        >>> hc = hec.HECClient()  # doctest: +REMOTE_DATA
        >>> print(hc.get_table_names())   # doctest: +REMOTE_DATA
        [(b'timed_see_flare',) (b'hi_event',) (b'yohkoh_flare_list',)
         (b'wind_mfi_bs_crossing_time',) (b'seeds_soho',) (b'seeds_stb',)
         ...
         (b'rhessi_hxr_flare',) (b'cactus_soho_flow',) (b'cactus_soho_cme',)
         (b'stereob_het_sep',)]

        """
        results = self.hec_client.service.getTableNames()
        tables = votable_handler(etree.tostring(results))
        return tables.array

    def select_table(self):
        """
        Creates a list of table names and prompts the user for a choice

        This takes the table of table names from get_table_names(), creates a
        list of the names, sorts them, then presents the tables in a
        convenient menu for the user to choose from. It returns a string
        containing the name of the table that the user picked.

        Returns
        -------
        `bytes`
            contains the name of the table that the user picked.

        Examples
        --------
        >>> from sunpy.net.helio import hec
        >>> hc = hec.HECClient()  # doctest: +REMOTE_DATA
        >>> hc.select_table()  # doctest: +REMOTE_DATA +SKIP

        """
        tables = self.get_table_names()
        table_list = [t[0] for t in tables if len(t[0]) > 0]
        table_list.sort()
        for index, table in enumerate(table_list):
            print(f'{index + 1} - {table.decode()}')

        while True:
            user_input = input(f"\nPlease enter a table number between 1 and {len(table_list)} "
                               "('e' to exit): ")
            if user_input.lower() == "e" or user_input.lower() == "exit":
                return None
            if user_input.isdigit() and 1 <= int(user_input) <= len(table_list):
                table_no = int(user_input)
                return table_list[table_no - 1]
            else:
                print(f"Input must be an integer between 1 and {len(table_list)}")
