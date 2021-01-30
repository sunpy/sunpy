"""
This module provides the `Fido
<sunpy.net.fido_factory.UnifiedDownloaderFactory>` instance of
`sunpy.net.fido_factory.UnifiedDownloaderFactory` it also provides the
`~sunpy.net.fido_factory.UnifiedResponse` class which
`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>` returns and the
`~sunpy.net.fido_factory.DownloadResponse` class that is returned by
`Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.

"""
import os
from pathlib import Path
from textwrap import dedent
from collections.abc import Sequence

import parfive

from astropy.table import Table

from sunpy.net import attr, vso
from sunpy.net.base_client import BaseClient, QueryResponseColumn, QueryResponseRow, QueryResponseTable
from sunpy.util.datatype_factory_base import BasicRegistrationFactory, NoMatchError
from sunpy.util.decorators import deprecated
from sunpy.util.parfive_helpers import Downloader, Results
from sunpy.util.util import get_width

__all__ = ['Fido', 'UnifiedResponse', 'UnifiedDownloaderFactory']


class UnifiedResponse(Sequence):
    """
    The object used to store results from `~sunpy.net.UnifiedDownloaderFactory.search`.

    The `~sunpy.net.Fido` object returns results from multiple different
    clients. So it is always possible to sub-select these results, you can
    index this object with two indices. The first index is the client index,
    i.e. corresponding to the results from the `~sunpy.net.vso.VSOClient`. The
    second index can be used to select records from the results returned from
    that client, for instance if you only want every second result you could
    index the second dimension with ``::2``.
    """

    def __init__(self, *results):
        """
        Parameters
        ----------
        results : `sunpy.net.base_client.QueryResponseTable`
            One or more QueryResponse objects.
        """
        self._list = []
        self._numfile = 0
        for result in results:
            if isinstance(result, QueryResponseRow):
                result = result.as_table()

            if isinstance(result, QueryResponseColumn):
                result = result.as_table()

            if not isinstance(result, QueryResponseTable):
                raise TypeError(
                    f"{type(result)} is not derived from sunpy.net.base_client.QueryResponseTable")

            self._list.append(result)
            self._numfile += len(result)

    def __len__(self):
        return len(self._list)

    def _getitem_string(self, aslice):
        ret = []
        for res in self._list:
            clientname = res.client.__class__.__name__
            if aslice.lower() == clientname.lower().split('client')[0]:
                ret.append(res)

        if len(ret) == 1:
            ret = ret[0]
        elif len(ret) == 0:
            raise IndexError(f"{aslice} is not a valid key, valid keys are: {','.join(self.keys())}")

        return ret

    def __getitem__(self, aslice):
        """
        Support slicing the UnifiedResponse as a 2D object.

        The first index is to the client and the second index is the records
        returned from those clients.
        """
        if isinstance(aslice, (int, slice)):
            ret = self._list[aslice]

        # using the client's name for indexing the responses.
        elif isinstance(aslice, str):
            ret = self._getitem_string(aslice)

        # Make sure we only have a length two slice.
        elif isinstance(aslice, tuple):
            if len(aslice) > 2:
                raise IndexError("UnifiedResponse objects can only "
                                 "be sliced with one or two indices.")

            # Indexing both client and records, but only for one client.
            if isinstance(aslice[0], str):
                intermediate = self._getitem_string(aslice[0])
            else:
                intermediate = self._list[aslice[0]]

            if isinstance(intermediate, list):
                ret = []
                for client_resp in intermediate:
                    ret.append(client_resp[aslice[1]])
            else:
                ret = intermediate[aslice[1]]

        else:
            raise IndexError("UnifiedResponse objects must be sliced with integers or strings.")

        if isinstance(ret, (QueryResponseTable, QueryResponseColumn, QueryResponseRow)):
            return ret

        return UnifiedResponse(*ret)

    @deprecated("2.1", "The same behaviour can now be obtained by indexing the object directly")
    def get_response(self, i):  # pragma: no cover
        """
        Get the actual response rather than another UnifiedResponse object.
        """
        return self._list[i]

    @deprecated("2.1", alternative="path_format_keys")
    def response_block_properties(self):
        return self.path_format_keys()

    def path_format_keys(self):
        """
        Returns all the names that can be used to format filenames.

        Only the keys which can be used to format all results from all
        responses contained in this `~.UnifiedResponse` are returned. Each
        individual response might have more keys available.

        Each one corresponds to a single column in the table, and the format
        syntax should match the dtype of that column, i.e. for a ``Time``
        object or a ``Quantity``.
        """
        s = self[0].path_format_keys()
        for table in self[1:]:
            s = s.intersection(table.path_format_keys())
        return s

    @property
    @deprecated("2.1", "All objects contained in the sequence are now Table objects")
    def tables(self):
        """
        Returns a list of `astropy.table.Table` for all responses present in a specific
        `~sunpy.net.fido_factory.UnifiedResponse` object. They can then be used
        to perform key-based indexing of objects of either type
        `sunpy.net.dataretriever.client.QueryResponse`, `sunpy.net.vso.QueryResponse` or
        `sunpy.net.jsoc.JSOCClient`

        Returns
        -------
        `list`
            A list of `astropy.table.Table`, consisting of data either from the
            `sunpy.net.dataretriever.client.QueryResponse`, `sunpy.net.vso.QueryResponse` or
            `sunpy.net.jsoc.JSOCClient`.
        """
        return list(self)

    @property
    @deprecated("2.1", "The same behaviour can be obtained by iterating over the object directly")
    def responses(self):  # pragma: no cover
        """
        A generator of all the `sunpy.net.dataretriever.client.QueryResponse`
        objects contained in the `~sunpy.net.fido_factory.UnifiedResponse`
        object.
        """
        for table in self:
            yield table

    def keys(self):
        """
        Names of the contained responses.

        One name may map to more than one response.
        """
        ret = []
        for res in self._list:
            clientname = res.client.__class__.__name__.lower().split('client')[0]
            if clientname not in ret:
                ret.append(clientname)
        return ret

    @property
    def file_num(self):
        """
        The number of records returned in all responses.
        """
        return self._numfile

    def _repr_html_(self):
        nprov = len(self)
        if nprov == 1:
            ret = 'Results from {} Provider:</br></br>'.format(len(self))
        else:
            ret = 'Results from {} Providers:</br></br>'.format(len(self))
        for block in self:
            ret += "{} Results from the {}:</br>".format(len(block),
                                                         block.client.__class__.__name__)
            ret += block._repr_html_()
            ret += '</br>'

        return ret

    def __repr__(self):
        return object.__repr__(self) + "\n" + str(self)

    def __str__(self):
        nprov = len(self)
        if nprov == 1:
            ret = 'Results from {} Provider:\n\n'.format(len(self))
        else:
            ret = 'Results from {} Providers:\n\n'.format(len(self))
        for block in self:
            ret += "{} Results from the {}:\n".format(len(block), block.client.__class__.__name__)
            lines = repr(block).split('\n')
            ret += '\n'.join(lines[1:])
            ret += '\n\n'

        return ret

    def show(self, *cols):
        """
        Displays response tables with desired columns for the Query.

        Parameters
        ----------
        \\*cols : `tuple`
            Name of columns to be shown.

        Returns
        -------
        `list` of `astropy.table.Table`
            A list of tables showing values for specified columns.
        """
        return type(self)(*[i.show(*cols) for i in self._list])

    @property
    def all_colnames(self):
        """
        Returns all the colnames in any of the tables in this response.

        Any column names in this list are valid inputs to :meth:`.UnifiedResponse.show`.
        """
        colnames = set(self[0].colnames)
        for resp in self[1:]:
            colnames.union(resp.colnames)
        return sorted(list(colnames))


query_walker = attr.AttrWalker()
"""
We construct an `AttrWalker` which calls `_make_query_to_client` for each
logical component of the query, i.e. any block which are ANDed together.
"""


@query_walker.add_creator(attr.DataAttr)
def _create_data(walker, query, factory):
    return factory._make_query_to_client(query)


@query_walker.add_creator(attr.AttrAnd)
def _create_and(walker, query, factory):
    return factory._make_query_to_client(*query.attrs)


@query_walker.add_creator(attr.AttrOr)
def _create_or(walker, query, factory):
    qblocks = []
    for attrblock in query.attrs:
        qblocks += walker.create(attrblock, factory)

    return qblocks


class UnifiedDownloaderFactory(BasicRegistrationFactory):
    """
    Fido is a unified data search and retrieval tool.

    It provides simultaneous access to a variety of online data sources, some
    cover multiple instruments and data products like the Virtual Solar
    Observatory and some are specific to a single source.

    For details of using `~sunpy.net.Fido` see :ref:`fido_guide`.

    """

    def search(self, *query):
        """
        Query for data in form of multiple parameters.

        Examples
        --------
        Query for LYRA timeseries data for the time range ('2012/3/4','2012/3/6')

        >>> from sunpy.net import Fido, attrs as a
        >>> import astropy.units as u
        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra) # doctest: +REMOTE_DATA

        Query for data from Nobeyama Radioheliograph and RHESSI

        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'),
        ...     (a.Instrument.norh & a.Wavelength(17*u.GHz)) | a.Instrument.rhessi)  # doctest: +REMOTE_DATA

        Query for 304 Angstrom SDO AIA data with a cadence of 10 minutes

        >>> import astropy.units as u
        >>> from sunpy.net import Fido, attrs as a
        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'),
        ...                        a.Instrument.aia,
        ...                        a.Wavelength(304*u.angstrom, 304*u.angstrom),
        ...                        a.Sample(10*u.minute))  # doctest: +REMOTE_DATA

        Parameters
        ----------
        query : `sunpy.net.vso.attrs`, `sunpy.net.jsoc.attrs`
            A query consisting of multiple parameters which define the
            requested data.  The query is specified using attributes from the
            VSO and the JSOC.  The query can mix attributes from the VSO and
            the JSOC.

        Returns
        -------
        `sunpy.net.fido_factory.UnifiedResponse`
            Container of responses returned by clients servicing query.

        Notes
        -----
        The conjunction 'and' transforms query into disjunctive normal form
        ie. query is now of form A & B or ((A & B) | (C & D))
        This helps in modularising query into parts and handling each of the
        parts individually.
        """
        query = attr.and_(*query)
        results = query_walker.create(query, self)

        # If we have searched the VSO but no results were returned, but another
        # client generated results, we drop the empty VSO results for tidiness.
        # This is because the VSO _can_handle_query is very broad because we
        # don't know the full list of supported values we can search for (yet).
        if len(results) > 1:
            vso_results = list(filter(lambda r: isinstance(r, vso.VSOQueryResponseTable), results))
            for vres in vso_results:
                if len(vres) == 0:
                    results.remove(vres)

        return UnifiedResponse(*results)

    def fetch(self, *query_results, path=None, max_conn=5, progress=True,
              overwrite=False, downloader=None, **kwargs):
        """
        Download the records represented by `~sunpy.net.base_client.QueryResponseTable` or
        `~sunpy.net.fido_factory.UnifiedResponse` objects.

        Parameters
        ----------
        query_results : `sunpy.net.fido_factory.UnifiedResponse` or `~sunpy.net.base_client.QueryResponseTable`
            Container returned by query method, or multiple.
        path : `str`
            The directory to retrieve the files into. Can refer to any fields
            in `UnifiedResponse.response_block_properties` via string formatting,
            moreover the file-name of the file downloaded can be referred to as file,
            e.g. "{source}/{instrument}/{time.start}/{file}".
        max_conn : `int`, optional
            The number of parallel download slots.
        progress : `bool`, optional
            If `True` show a progress bar showing how many of the total files
            have been downloaded. If `False`, no progress bars will be shown at all.
        overwrite : `bool` or `str`, optional
            Determine how to handle downloading if a file already exists with the
            same name. If `False` the file download will be skipped and the path
            returned to the existing file, if `True` the file will be downloaded
            and the existing file will be overwritten, if ``'unique'`` the filename
            will be modified to be unique.
        downloader : `parfive.Downloader`, optional
            The download manager to use. If specified the ``max_conn``,
            ``progress`` and ``overwrite`` arguments are ignored.

        Returns
        -------
        `parfive.Results`

        Examples
        --------
        >>> from sunpy.net.attrs import Time, Instrument
        >>> unifresp = Fido.search(Time('2012/3/4','2012/3/5'), Instrument('EIT'))  # doctest: +REMOTE_DATA
        >>> filepaths = Fido.fetch(unifresp)  # doctest: +SKIP

        If any downloads fail, they can be retried by passing the `parfive.Results` object back into ``fetch``.

        >>> filepaths = Fido.fetch(filepaths)  # doctest: +SKIP

        """
        if path is not None:
            exists = list(filter(lambda p: p.exists(), Path(path).resolve().parents))

            if not os.access(exists[0], os.W_OK):
                raise PermissionError('You do not have permission to write'
                                      f' to the directory {exists[0]}.')

        if "wait" in kwargs:
            raise ValueError("wait is not a valid keyword argument to Fido.fetch.")

        if downloader is None:
            downloader = Downloader(max_conn=max_conn, progress=progress, overwrite=overwrite)
        elif not isinstance(downloader, parfive.Downloader):
            raise TypeError("The downloader argument must be a parfive.Downloader object.")

        # Handle retrying failed downloads
        retries = [isinstance(arg, Results) for arg in query_results]
        if all(retries):
            results = Results()
            for retry in query_results:
                dr = downloader.retry(retry)
                results.data += dr.data
                results._errors += dr._errors
            return results
        elif any(retries):
            raise TypeError("If any arguments to fetch are "
                            "`parfive.Results` objects, all arguments must be.")

        reslist = []
        for query_result in query_results:
            if isinstance(query_result, QueryResponseRow):
                responses = [query_result.as_table()]
            elif isinstance(query_result, QueryResponseTable):
                responses = [query_result]
            elif isinstance(query_result, UnifiedResponse):
                responses = query_result
            else:
                raise ValueError(f"Query result has an unrecognized type: {type(query_result)} "
                                 "Allowed types are QueryResponseRow, QueryResponseTable or UnifiedResponse.")
            for block in responses:
                result = block.client.fetch(block, path=path,
                                            downloader=downloader,
                                            wait=False, **kwargs)
                if result is not NotImplemented:
                    reslist.append(result)

        results = downloader.download()
        # Combine the results objects from all the clients into one Results
        # object.
        for result in reslist:
            if result is None:
                continue
            if not isinstance(result, Results):
                raise TypeError(
                    "If wait is False a client must return a parfive.Downloader and either None"
                    " or a parfive.Results object.")
            results.data += result.data
            results._errors += result.errors

        return results

    def __call__(self, *args, **kwargs):
        raise TypeError(f"'{self.__class__.__name__}' object is not callable")

    def _check_registered_widgets(self, *args):
        """Factory helper function"""
        candidate_widget_types = list()
        for key in self.registry:
            if self.registry[key](*args):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)
        if n_matches == 0:
            # There is no default client
            raise NoMatchError("This query was not understood by any clients. Did you miss an OR?")

        return candidate_widget_types

    def _make_query_to_client(self, *query):
        """
        Given a query, look up the client and perform the query.

        Parameters
        ----------
        query : collection of `~sunpy.net.vso.attr` objects

        Returns
        -------
        results : `list`

        client : `object`
            Instance of client class
        """
        candidate_widget_types = self._check_registered_widgets(*query)
        results = []
        for client in candidate_widget_types:
            tmpclient = client()
            kwargs = dict()
            # Handle the change in response format in the VSO
            if isinstance(tmpclient, vso.VSOClient):
                kwargs = dict(response_format="table")
            results.append(tmpclient.search(*query, **kwargs))

        # This method is called by `search` and the results are fed into a
        # UnifiedResponse object.
        return results

    def __repr__(self):
        return object.__repr__(self) + "\n" + str(self)

    def __str__(self):
        """
        This enables the "pretty" printing of the Fido Clients.
        """
        return self._print_clients()

    def _repr_html_(self):
        """
        This enables the "pretty" printing of the Fido Clients with html.
        """
        return self._print_clients(html=True)

    def _print_clients(self, html=False) -> str:
        width = -1 if html else get_width()

        t = Table(names=["Client", "Description"], dtype=["U80", "U120"])
        lines = ["sunpy.net.Fido", dedent(self.__doc__)]
        if html:
            lines = [f"<p>{line}</p>" for line in lines]
        for key in BaseClient._registry.keys():
            t.add_row((key.__name__, dedent(
                key.__doc__.partition("\n\n")[0].replace("\n    ", " "))))
        lines.extend(t.pformat_all(show_dtype=False, max_width=width, align="<", html=html))
        return '\n'.join(lines)


Fido = UnifiedDownloaderFactory(
    registry=BaseClient._registry, additional_validation_functions=['_can_handle_query'])
