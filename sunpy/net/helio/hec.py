"""
Access the Helio Event Catalogue
"""
import io
import os

from lxml import etree
from requests import Session
from requests.exceptions import SSLError
from zeep import Client
from zeep.transports import Transport

from astropy.io.votable.table import parse_single_table

from sunpy import log
from sunpy.net import attrs as a
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.helio import attrs as ha
from sunpy.net.helio import parser
from sunpy.time import parse_time
from sunpy.util.exceptions import warn_user

__all__ = ['HECClient', 'HECResponse']


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
    for i in range(len(votable.array)):
        item = votable.array[i][0]
        if isinstance(item, bytes):
            votable.array[i] = (votable.array[i][0].decode(),)
    fake_file.close()
    return votable


class HECResponse(QueryResponseTable):
    """
    A container for data returned from HEC searches.
    """


class HECClient(BaseClient):
    """
    Provides access to the HELIO webservices.
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
        >>> hc = hec.HECClient()  # doctest: +SKIP
        """
        if link is None:
            # The default wsdl file
            link = parser.wsdl_retriever()
        session = Session()
        # This is for use in our test suite.
        session.verify = not (bool(os.environ.get("NO_VERIFY_HELIO_SSL", 0)))
        transport = Transport(session=session)
        try:
            self.hec_client = Client(link, transport=transport)
        except SSLError:
            log.warning('SSL verification error for HEC client.\n'
                     'Set the \'NO_VERIFY_HELIO_SSL\' environment variable disable SSL verification for Helio.')
            self.hec_client = None

    @classmethod
    def _can_handle_query(cls, *query):
        required = {a.Time}
        optional = {ha.MaxRecords, ha.TableName}
        return cls.check_attr_types_in_query(query, required, optional)

    @classmethod
    def _attrs_module(cls):
        return 'helio', 'sunpy.net.helio.attrs'

    def search(self, *args, **kwargs):
        """
        The simple interface to query the wsdl service.

        Used to utilize the service's TimeQuery() method, this is a simple
        interface between the sunpy module library and the web-service's API.

        .. note::
           By default the maximum records returned by the service are limited to 500.
           To obtain more results ``a.helio.MaxRecords`` must be set to a higher value.

        Examples
        --------
        >>> from sunpy.net.helio import attrs as ha
        >>> from sunpy.net import attrs as a, Fido
        >>> timerange = a.Time('2005/01/03', '2005/12/03')
        >>> res = Fido.search(timerange, ha.MaxRecords(10),
        ...                   ha.TableName('rhessi_hxr_flare'))  # doctest: +SKIP
        >>> res   # doctest: +SKIP
        <sunpy.net.fido_factory.UnifiedResponse object at ...>
        Results from 1 Provider:
        <BLANKLINE>
        10 Results from the HECClient:
        hec_id      time_start          time_peak      ... energy_kev flare_number
        ------ ------------------- ------------------- ... ---------- ------------
         31463 2005-01-03T01:37:36 2005-01-03T01:37:54 ...          6      5010320
         31464 2005-01-03T01:51:36 2005-01-03T01:59:18 ...         12      5010301
         31465 2005-01-03T03:26:28 2005-01-03T03:42:50 ...          6      5010332
         31466 2005-01-03T03:46:04 2005-01-03T04:07:10 ...         12      5010302
         31467 2005-01-03T05:00:24 2005-01-03T05:00:30 ...          6      5010313
         31468 2005-01-03T06:40:48 2005-01-03T06:42:46 ...          6      5010314
         31469 2005-01-03T08:27:56 2005-01-03T08:28:26 ...          6      5010334
         31470 2005-01-03T09:31:00 2005-01-03T09:33:34 ...          6      5010322
         31471 2005-01-03T09:34:52 2005-01-03T09:59:46 ...          6      5010336
         31472 2005-01-03T11:06:48 2005-01-03T11:07:18 ...         12      5010304
        <BLANKLINE>
        <BLANKLINE>
        """
        if self.hec_client is None:
            table = HECResponse([], client=self)
            return table
        qrdict = {}
        for elem in args:
            if isinstance(elem, a.Time):
                qrdict['Time'] = elem
            elif isinstance(elem, ha.MaxRecords):
                qrdict['max_records'] = elem.value
            elif isinstance(elem, ha.TableName):
                qrdict['table_name'] = elem.value
            else:
                raise ValueError(
                    f"{elem.__class__.__name__} should be a ``attrs.Time``, ``attrs.hek.MaxRecords`` or ``attrs.hek.TableName`` attribute.")
        qrdict.update(kwargs)
        table = qrdict.get('table_name', None)
        start_time = qrdict['Time'].start
        end_time = qrdict['Time'].end
        max_records = qrdict.get('max_records', 500)
        while table is None:
            table = self.select_table()
        start_time = parse_time(start_time)
        end_time = parse_time(end_time)
        results = self.hec_client.service.TimeQuery(STARTTIME=start_time.isot,
                                                    ENDTIME=end_time.isot,
                                                    FROM=table,
                                                    MAXRECORDS=max_records)
        results = votable_handler(etree.tostring(results))
        table = HECResponse(results.to_table(), client=self)
        if len(table) == max_records:
            warn_user("Number of results is the same as current limit. "
                      "It is possible your query has been truncated. "
                      "If you want to change this, set `a.helio.MaxRecords` to a higher value.")
        return table

    def get_table_names(self):
        """
        Returns a list of the available tables to query.

        Returns the names of all the tables that can be queried via the
        webservice.

        Returns
        -------
        tables.array: `numpy.ma.core.MaskedArray`
            A VOtable table of available tables names.

        Examples
        --------
        >>> from sunpy.net.helio import hec
        >>> hc = hec.HECClient()  # doctest: +SKIP
        >>> print(hc.get_table_names())  # doctest: +SKIP
        [('timed_see_flare',) ('hi_event',) ('yohkoh_flare_list',)
         ('wind_mfi_bs_crossing_time',) ('seeds_soho',) ('seeds_stb',)
         ...
         ('rhessi_hxr_flare',) ('cactus_soho_flow',) ('cactus_soho_cme',)
         ('stereob_het_sep',)]
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
        `str`
            Contains the name of the table that the user picked.

        Examples
        --------
        >>> from sunpy.net.helio import hec  # doctest: +SKIP
        >>> hc = hec.HECClient()  # doctest: +SKIP
        >>> hc.select_table()  # doctest: +SKIP
        """
        tables = self.get_table_names()
        table_list = sorted([t[0] for t in tables if len(t[0]) > 0])
        for index, table in enumerate(table_list):
            print(f'{index + 1} - {table}')
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

    def fetch(self, *args, **kwargs):
        """
        This is a no operation function as this client does not download data.
        """
        return NotImplemented
