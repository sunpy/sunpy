"""
This module contains the factory code for the ``sunpy.net.Fido`` factory. As
well as defining ``Fido`` which is an instance of
`sunpy.net.fido_factory.UnifiedDownloaderFactory` it defines
`~sunpy.net.fido_factory.UnifiedResponse` which ``Fido.search`` returns.
`~sunpy.net.fido_factory.DownloadResponse` objects are returned by ``Fido.fetch``.

"""
# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014
from collections import MutableSequence

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError

from sunpy.net.dataretriever.clients import CLIENTS
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.vso import VSOClient
from . import attr
from . import attrs as a

__all__ = ['Fido', 'UnifiedResponse', 'UnifiedDownloaderFactory', 'DownloadResponse']


class UnifiedResponse(MutableSequence):
    """
    The object used to store responses from the unified downloader.
    """
    def __init__(self, lst):
        """
        Input to this constructor can be one of a few things:

        1. A list of one UnifiedResponse object
        2. A list of tuples (QueryResponse, client)
        """

        tmplst = []
        # numfile is the number of files not the number of results.
        self._numfile = 0
        if isinstance(lst, QueryResponse):
            if not hasattr(lst, 'client'):
                raise("QueryResponse is only a valid input if it has a client attribute.")
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
                    raise Exception("{} is not a valid input to UnifiedResponse.".format(type(lst)))

        self._list = tmplst

    def __len__(self):
        return len(self._list)

    def __getitem__(self, aslice):
        ret = self._list[aslice]
        if ret:
            if isinstance(ret, list):
                return type(self)(ret)
            else:
                return type(self)(ret)

        return ret

    def __delitem__(self, i):
        del self._list[i]

    def __setitem__(self, aslice, v):
        self._list[aslice] = v

    def __iter__(self):
        return self.responses

    def insert(self, i, v):
        self._list.insert(i, v)

    def get_response(self, i):
        """
        Get the actual response rather than another UnifiedResponse object.
        """
        return self._list[i]

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
        ret = ''
        for block in self.responses:
            ret += block._repr_html_()

        return ret


class DownloadResponse(list):
    """
    Object returned by clients servicing the query.
    """

    def __init__(self, lst):
        super(DownloadResponse, self).__init__(lst)

    def wait(self, progress=True):
        """
        Waits for all files to download completely and then return.

        Parameters
        ----------
        progress : `bool`
            if true, display a progress bar.

        Returns
        -------
        List of file paths to which files have been downloaded.
        """
        filelist = []
        for resobj in self:
            filelist.extend(resobj.wait(progress=progress))

        return filelist


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
    att = {type(x) for x in query.attrs}
    if a.Time not in att:
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
    sunpy.net.Fido(\*args, \*\*kwargs)

    Search and Download data from a variety of supported sources.
    """
    def search(self, *query):
        """
        Query for data in form of multiple parameters.

        Examples
        --------
        Query for LYRALightCurve data for the time range ('2012/3/4','2012/3/6')

        >>> from sunpy.net.vso.attrs import Time, Instrument
        >>> unifresp = Fido.search(Time('2012/3/4', '2012/3/6'), Instrument('lyra'))

        Query for data from Nobeyama Radioheliograph and RHESSI

        >>> unifresp = Fido.search(Time('2012/3/4', '2012/3/6'), Instrument('norh') | Instrument('rhessi'))

        Query for 304 Angstrom SDO AIA data with a cadence of 10 minutes

        >>> import astropy.units as u
        >>> from sunpy.net.vso.attrs import Time, Instrument, Wavelength, Sample
        >>> unifresp = Fido.search(Time('2012/3/4', '2012/3/6'), Instrument('AIA'), Wavelength(304*u.angstrom, 304*u.angstrom), Sample(10*u.minute))

        Parameters
        ----------
        query : `sunpy.net.vso.attrs`, `sunpy.net.jsoc.attrs`
            A query consisting of multiple parameters which define the
            requested data.  The query is specified using attributes from the
            VSO and the JSOC.  The query can mix attributes from the VSO and
            the JSOC.

        Returns
        -------
        `sunpy.net.fido_factory.UnifiedResponse` object
            Container of responses returned by clients servicing query.

        Notes
        -----
        The conjunction 'and' transforms query into disjunctive normal form
        ie. query is now of form A & B or ((A & B) | (C & D))
        This helps in modularising query into parts and handling each of the
        parts individually.
        """
        query = attr.and_(*query)
        return UnifiedResponse(query_walker.create(query, self))

    def fetch(self, query_result, wait=True, progress=True, **kwargs):
        """
        Downloads the files pointed at by URLs contained within UnifiedResponse
        object.

        Parameters
        ----------
        query_result : `sunpy.net.fido_factory.UnifiedResponse`
            Container returned by query method.

        wait : `bool`
            fetch will wait until the download is complete before returning.

        progress : `bool`
            Show a progress bar while the download is running.

        Returns
        -------
        `sunpy.net.fido_factory.DownloadResponse`

        Example
        --------
        >>> from sunpy.net.vso.attrs import Time, Instrument
        >>> unifresp = Fido.search(Time('2012/3/4','2012/3/6'), Instrument('AIA'))
        >>> downresp = Fido.get(unifresp)
        >>> file_paths = downresp.wait()
        """
        reslist = []
        for block in query_result.responses:
            reslist.append(block.client.get(block, **kwargs))

        results = DownloadResponse(reslist)

        if wait:
            return results.wait(progress=progress)
        else:
            return results

    def __call__(self, *args, **kwargs):
        raise TypeError("'{}' object is not callable".format(
            self.__class__.__name__))

    def _check_registered_widgets(self, *args):
        """Factory helper function"""
        candidate_widget_types = list()
        for key in self.registry:

            if self.registry[key](*args):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)
        if n_matches == 0:
            # There is no default client
            raise NoMatchError(
                "This query was not understood by any clients. Did you miss an OR?")
        elif n_matches == 2:
            # If two clients have reported they understand this query, and one
            # of them is the VSOClient, then we ignore VSOClient.
            if VSOClient in candidate_widget_types:
                candidate_widget_types.remove(VSOClient)

        # Finally check that we only have one match.
        if len(candidate_widget_types) > 1:
            candidate_names = [cls.__name__ for cls in candidate_widget_types]
            raise MultipleMatchError(
                "The following clients matched this query. "
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

        client : Instance of client class
        """
        candidate_widget_types = self._check_registered_widgets(*query)
        tmpclient = candidate_widget_types[0]()
        return tmpclient.query(*query), tmpclient


Fido = UnifiedDownloaderFactory(registry=CLIENTS,
                                additional_validation_functions=['_can_handle_query'])
