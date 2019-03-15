"""
This module provides the `Fido
<sunpy.net.fido_factory.UnifiedDownloaderFactory>` instance of
`sunpy.net.fido_factory.UnifiedDownloaderFactory` it also provides the
`~sunpy.net.fido_factory.UnifiedResponse` class which
`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>` returns and the
`~sunpy.net.fido_factory.DownloadResponse` class that is returned by
`Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.

"""
from collections import Sequence

from parfive import Downloader, Results

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError

from sunpy.net.base_client import BaseClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.vso import VSOClient, QueryResponse as vsoQueryResponse

from sunpy.net import attr
from sunpy.net import attrs as a

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

    def __init__(self, lst):
        """
        Parameters
        ----------
        lst : `object`
            A single instance or an iterable of ``(QueryResponse, client)``
            pairs or ``QueryResponse`` objects with a ``.client`` attribute.
        """

        tmplst = []
        # numfile is the number of files not the number of results.
        self._numfile = 0
        if isinstance(lst, (QueryResponse, vsoQueryResponse)):
            if not hasattr(lst, 'client'):
                raise ValueError(
                    ("A {} object is only a valid input to UnifiedResponse "
                     "if it has a client attribute.").
                    format(type(lst).__name__))
            tmplst.append(lst)
            self._numfile = len(lst)
        else:
            for block in lst:
                if isinstance(block, tuple) and len(block) == 2:
                    block[0].client = block[1]
                    tmplst.append(block[0])
                    self._numfile += len(block[0])
                elif hasattr(block, 'client'):
                    tmplst.append(block)
                    self._numfile += len(block)
                else:
                    raise ValueError(
                        "{} is not a valid input to UnifiedResponse.".format(type(lst)))
        self._list = tmplst

    def __len__(self):
        return len(self._list)

    def __iter__(self):
        return self.responses

    def _handle_record_slice(self, client_resp, record_slice):
        """
        Given a slice to be applied to the results from a single client, return
        an object of the same type as client_resp.
        """
        # When we subindex, we want to persist the type of the response object.
        resp_type = type(client_resp)

        # Make sure we always have an iterable, as most of the response objects
        # expect one.
        if isinstance(record_slice, int):
            resp = [client_resp[record_slice]]
        else:
            resp = client_resp[record_slice]

        # Reconstruct a response object with the sub-indexed records.
        ret = resp_type(resp)
        # Make sure we pass the client back out again.
        ret.client = client_resp.client

        return ret

    def __getitem__(self, aslice):
        """
        Support slicing the UnifiedResponse as a 2D object.

        The first index is to the client and the second index is the records
        returned from those clients.
        """
        # Just a single int as a slice, we are just indexing client.
        if isinstance(aslice, (int, slice)):
            ret = self._list[aslice]

        # Make sure we only have a length two slice.
        elif isinstance(aslice, tuple):
            if len(aslice) > 2:
                raise IndexError("UnifiedResponse objects can only "
                                 "be sliced with one or two indices.")

            # Indexing both client and records, but only for one client.
            if isinstance(aslice[0], int):
                client_resp = self._list[aslice[0]]
                ret = self._handle_record_slice(client_resp, aslice[1])

            # Indexing both client and records for multiple clients.
            else:
                intermediate = self._list[aslice[0]]
                ret = []
                for client_resp in intermediate:
                    resp = self._handle_record_slice(client_resp, aslice[1])
                    ret.append(resp)

        else:
            raise IndexError("UnifiedResponse objects must be sliced with integers.")

        return UnifiedResponse(ret)

    def get_response(self, i):
        """
        Get the actual response rather than another UnifiedResponse object.
        """
        return self._list[i]

    def response_block_properties(self):
        """
        Returns a set of class attributes on all the response blocks.

        Returns
        -------
        s : list
            List of strings, containing attribute names in the response blocks.
        """
        s = self.get_response(0).response_block_properties()
        for i in range(1, len(self)):
            s.intersection(self.get_response(i).response_block_properties())
        return s

    @property
    def responses(self):
        """
        A generator of all the `sunpy.net.dataretriever.client.QueryResponse`
        objects contained in the `~sunpy.net.fido_factory.UnifiedResponse`
        object.
        """
        for i in range(len(self)):
            yield self.get_response(i)

    @property
    def file_num(self):
        return self._numfile

    def _repr_html_(self):
        nprov = len(self)
        if nprov == 1:
            ret = 'Results from {} Provider:</br></br>'.format(len(self))
        else:
            ret = 'Results from {} Providers:</br></br>'.format(len(self))
        for block in self.responses:
            ret += "{} Results from the {}:</br>".format(len(block),
                                                         block.client.__class__.__name__)
            ret += block._repr_html_()
            ret += '</br>'

        return ret

    def __repr__(self):
        ret = super(UnifiedResponse, self).__repr__()
        ret += '\n' + str(self)

        return ret

    def __str__(self):
        nprov = len(self)
        if nprov == 1:
            ret = 'Results from {} Provider:\n\n'.format(len(self))
        else:
            ret = 'Results from {} Providers:\n\n'.format(len(self))
        for block in self.responses:
            ret += "{} Results from the {}:\n".format(len(block), block.client.__class__.__name__)
            lines = repr(block).split('\n')
            ret += '\n'.join(lines[1:])
            ret += '\n\n'

        return ret


"""
Construct a simple AttrWalker to split up searches into blocks of attrs being
'anded' with AttrAnd.

This pipeline only understands AttrAnd and AttrOr, Fido.search passes in an
AttrAnd object of all the query parameters, if an AttrOr is encountered the
query is split into the component parts of the OR, which at somepoint will end
up being an AttrAnd object, at which point it is passed into
_get_registered_widget.
"""
query_walker = attr.AttrWalker()


@query_walker.add_creator(attr.AttrAnd)
def _create_and(walker, query, factory):
    is_time = any([isinstance(x, a.Time) for x in query.attrs])
    if not is_time:
        error = "The following part of the query did not have a time specified:\n"
        for at in query.attrs:
            error += str(at) + ', '
        raise ValueError(error)

    # Return the response and the client
    return [factory._make_query_to_client(*query.attrs)]


@query_walker.add_creator(attr.AttrOr)
def _create_or(walker, query, factory):
    qblocks = []
    for attrblock in query.attrs:
        qblocks.extend(walker.create(attr.and_(attrblock), factory))

    return qblocks


class UnifiedDownloaderFactory(BasicRegistrationFactory):
    """
    sunpy.net.Fido(\\*args, \\*\\*kwargs)

    Search and Download data from a variety of supported sources.
    """

    def search(self, *query):
        """
        Query for data in form of multiple parameters.

        Examples
        --------
        Query for LYRALightCurve data for the time range ('2012/3/4','2012/3/6')

        >>> from sunpy.net import Fido, attrs as a
        >>> import astropy.units as u
        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('lyra')) # doctest: +REMOTE_DATA

        Query for data from Nobeyama Radioheliograph and RHESSI

        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'),
        ...     (a.Instrument('norh') & a.Wavelength(17*u.GHz)) | a.Instrument('rhessi'))  # doctest: +REMOTE_DATA

        Query for 304 Angstrom SDO AIA data with a cadence of 10 minutes

        >>> import astropy.units as u
        >>> from sunpy.net import Fido, attrs as a
        >>> unifresp = Fido.search(a.Time('2012/3/4', '2012/3/6'),
        ...                        a.Instrument('AIA'),
        ...                        a.Wavelength(304*u.angstrom, 304*u.angstrom),
        ...                        a.vso.Sample(10*u.minute))  # doctest: +REMOTE_DATA

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
        """  # noqa
        query = attr.and_(*query)
        return UnifiedResponse(query_walker.create(query, self))

    def fetch(self, *query_results, path=None, max_conn=5, progress=True,
              overwrite=False, downloader=None, **kwargs):
        """
        Download the records represented by
        `~sunpy.net.fido_factory.UnifiedResponse` objects.

        Parameters
        ----------
        query_results : `sunpy.net.fido_factory.UnifiedResponse`
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
            and the existing file will be overwritten, if `'unique'` the filename
            will be modified to be unique.

        downloader : `parfive.Downloader`, optional
            The download manager to use. If specified the ``max_conn``,
            ``progress`` and ``overwrite`` arguments are ignored.

        Returns
        -------
        `parfive.Results`

        Examples
        --------
        >>> from sunpy.net.vso.attrs import Time, Instrument
        >>> unifresp = Fido.search(Time('2012/3/4','2012/3/5'), Instrument('EIT'))  # doctest: +REMOTE_DATA
        >>> filepaths = Fido.fetch(unifresp)  # doctest: +SKIP

        If any downloads fail, they can be retried by passing the `parfive.Results` object back into ``fetch``.

        >>> filepaths = Fido.fetch(filepaths)  # doctest: +SKIP

        """

        if "wait" in kwargs:
            raise ValueError("wait is not a valid keyword argument to Fido.fetch.")

        if downloader is None:
            downloader = Downloader(max_conn=max_conn, progress=progress, overwrite=overwrite)
        elif not isinstance(downloader, Downloader):
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
            for block in query_result.responses:
                reslist.append(block.client.fetch(block, path=path,
                                                  downloader=downloader,
                                                  wait=False, **kwargs))

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
        raise TypeError("'{}' object is not callable".format(self.__class__.__name__))

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
        elif n_matches == 2:
            # If two clients have reported they understand this query, and one
            # of them is the VSOClient, then we ignore VSOClient.
            if VSOClient in candidate_widget_types:
                candidate_widget_types.remove(VSOClient)

        # Finally check that we only have one match.
        if len(candidate_widget_types) > 1:
            candidate_names = [cls.__name__ for cls in candidate_widget_types]
            raise MultipleMatchError("The following clients matched this query. "
                                     "Please make your query more specific.\n"
                                     "{}".format(candidate_names))

        return candidate_widget_types

    def _make_query_to_client(self, *query):
        """
        Given a query, look up the client and perform the query.

        Parameters
        ----------
        query : collection of `~sunpy.net.vso.attr` objects

        Returns
        -------
        response : `~sunpy.net.dataretriever.client.QueryResponse`

        client : `object`
            Instance of client class
        """
        candidate_widget_types = self._check_registered_widgets(*query)
        tmpclient = candidate_widget_types[0]()
        return tmpclient.search(*query), tmpclient


Fido = UnifiedDownloaderFactory(
    registry=BaseClient._registry, additional_validation_functions=['_can_handle_query'])
