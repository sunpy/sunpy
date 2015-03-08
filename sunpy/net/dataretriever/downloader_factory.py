#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError

from .. import attr
from .client import GenericClient

__all__ = ['Fido']

class UnifiedResponse(list):

    def __init__(self, lst):

        tmplst = []
        for block in lst:
            block[0].client = block[1]
            tmplst.append(block[0])
        super(UnifiedResponse, self).__init__(tmplst)
        self._numfile =0
        for qblock in self:
            self._numfile += len(qblock)

    @property
    def file_num(self):
        return self._numfile


class downloadresponse(list):
    """
    List of Results object returned by clients servicing the query.
    """
    def __init__(self,lst):

        super(downloadresponse, self).__init__(lst)

    def wait(self):
        """
        Waits for all files to download completely and then return.
        Returns
        -------
        list of file paths to which files have been downloaded.
        """
        filelist = []
        for resobj in self:
            filelist.extend(resobj.wait())

        return filelist


qwalker = attr.AttrWalker()

@qwalker.add_creator(attr.AttrAnd)
def _create(wlk, query, dobj):
    qresponseobj, qclient = dobj._get_registered_widget(*query.attrs)
    return [(qresponseobj, qclient)]


@qwalker.add_creator(attr.AttrOr)
def _create(wlk, query, dobj):
    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr, dobj))

    return qblocks


class UnifiedDownloaderFactory(BasicRegistrationFactory):

    def search(self, *query):
        """
        Query for data in form of multiple parameters.
        Examples
        --------
        Query for LYRALightCurve data from timerange('2012/3/4','2012/3/6')
        >>> unifresp = Fido.query(Time('2012/3/4','2012/3/6'),Instrument('lyra'))
        >>> unifresp = Fido.query(Time('2012/3/4','2012/3/6'),Instrument('norh') | Instrument('rhessi'))
        >>> unifresp = Fido.query(Time('2012/3/4','2012/3/6'),Instrument('AIA'),
                       Wave(304, 304),Sample(60*10))

        Parameters
        ----------
        query: Mutiple parameters,VSO-styled query. Attributes from JSOC, VSO both can be used.

        Returns
        -------
        UnifiedResponse object: Container of responses returned by clients servicing query.

        Notes
        -----
        and_ tranforms query into disjunctive normal form
        ie. query is now of form A & B or ((A & B) | (C & D))
        This helps in modularising query into parts and handling each of the parts individually.
        """
        query = attr.and_(*query)
        return UnifiedResponse(qwalker.create(query, self))

    def fetch(self, qr, **kwargs):
        """
        Downloads the files pointed at by URLS contained within UnifiedResponse Object.
        Parameters
        ----------
        qr : UnifiedResponse Object
            Container returned by query method.

        Returns
        -------
        DownloadResponse Object
            List of Results object with an additional wait method.

        Example
        --------
        >>> unifresp = Fido.query(Time('2012/3/4','2012/3/6'),Instrument('AIA'))
        >>> downresp = Fido.get(unifresp)
        >>> file_paths = downresp.wait()
        """
        reslist =[]
        for block in qr:
            reslist.append(block.client.get(block, **kwargs))

        return downloadresponse(reslist)

    def __call__(self, *args, **kwargs):
        pass


    def _check_registered_widgets(self, *args, **kwargs):
        """Factory helper function"""
        candidate_widget_types = list()
        for key in self.registry:

            if self.registry[key](*args):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)
        if n_matches == 0:
            if self.default_widget_type is None:
                raise NoMatchError("Query {0} can not be handled in its current form".format(args))
            else:
                return  [self.default_widget_type]
        elif n_matches > 1:
            # This is a hack, VSO services all Instruments.
            # TODO: VSOClient._can_handle_query should know what values of
            # Instrument VSO can handle.
            for candidate_client in candidate_widget_types:
                if issubclass(candidate_client, GenericClient):
                    return [candidate_client]

            candidate_names = [cls.__name__ for cls in candidate_widget_types]
            raise MultipleMatchError("Too many candidates clients can service your query {0}".format(candidate_names))

        return candidate_widget_types

    def _get_registered_widget(self, *args, **kwargs):
        """Factory helper function"""
        candidate_widget_types = self._check_registered_widgets(*args)
        tmpclient = candidate_widget_types[0]()
        return tmpclient.query(*args), tmpclient


Fido = UnifiedDownloaderFactory(additional_validation_functions = ['_can_handle_query'])
