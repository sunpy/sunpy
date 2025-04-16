# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014


from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient, QueryResponse

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']


class NOAAIndicesClient(GenericClient):
    """
    Provides access to the NOAA solar cycle indices.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed dataset so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_indices)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAIndicesClient:
    Source: https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json
    <BLANKLINE>
     Instrument     Physobs     Source Provider
    ------------ -------------- ------ --------
    NOAA-Indices sunspot number   SIDC     SWPC
    <BLANKLINE>
    <BLANKLINE>

    """
    required = {a.Instrument}

    @property
    def info_url(self):
        return 'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'

    def search(self, *args, **kwargs):
        rowdict = self._get_match_dict(*args, **kwargs)
        for key in rowdict:
            if isinstance(rowdict[key], list):
                # uses first value among the list of possible values corresponding to an Attr
                # returned by `get_match_dict()` to be shown in query response table.
                rowdict[key] = rowdict[key][0]
        rowdict['url'] = 'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'
        rowdict['Instrument'] = 'NOAA-Indices'
        # These results are not dependent on time, but we allow time as a
        # parameter for easy searching, so remove time from the results table
        # injected by GenericClient.
        rowdict.pop('Start Time', None)
        rowdict.pop('End Time', None)
        return QueryResponse([rowdict], client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Indices', 'Recent Solar Indices of Observed Monthly Mean Values')],
            attrs.Physobs: [('sunspot number', 'Sunspot Number.')],
            attrs.Source: [('SIDC', 'The Solar Influence Data Analysis Center')],
            attrs.Provider: [('SWPC', 'The Space Weather Prediction Center.')],
            attrs.Time: [('*')]}
        return adict


class NOAAPredictClient(GenericClient):
    """
    Provides access to the NOAA SWPC predicted sunspot Number and 10.7 cm radio flux values.

    Uses the `SWPC NOAA archive <https://services.swpc.noaa.gov/json/solar-cycle/>`__.
    This is a fixed prediction so the result is independent of the time range.

    Examples
    --------

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.noaa_predict)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the NOAAPredictClient:
    Source: https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json
    <BLANKLINE>
     Instrument     Physobs     Source Provider
    ------------ -------------- ------ --------
    NOAA-Predict sunspot number   ISES     SWPC
    <BLANKLINE>
    <BLANKLINE>

    """
    required = {a.Instrument}

    @property
    def info_url(self):
        return 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'

    def search(self, *args, **kwargs):
        rowdict = self._get_match_dict(*args, **kwargs)
        for key in rowdict:
            if isinstance(rowdict[key], list):
                # uses first value among the list of possible values corresponding to an Attr
                # returned by `get_match_dict()` to be shown in query response table.
                rowdict[key] = rowdict[key][0]
        rowdict['url'] = 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'
        rowdict['Instrument'] = 'NOAA-Predict'
        # These results are not dependent on time, but we allow time as a
        # parameter for easy searching, so remove time from the results table
        # injected by GenericClient.
        rowdict.pop('Start Time', None)
        rowdict.pop('End Time', None)
        return QueryResponse([rowdict], client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ('NOAA-Predict', 'Predicted Sunspot Number And Radio Flux Values With Expected Ranges.')],
            attrs.Physobs: [('sunspot number', 'Sunspot Number.')],
            attrs.Source: [('ISES', 'The International Space Environmental Services.')],
            attrs.Provider: [('SWPC', 'The Space Weather Prediction Center.')],
            attrs.Time: [('*')]}
        return adict


class SRSClient(GenericClient):
    """
    Provides access to the NOAA SWPC solar region summary data.

    `The data is retrieved from NOAA's NCEI archive
    <https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_region_summaries/>`__.

    Notes
    -----
    Data pre-1996 is in free-form text, which cannot be parsed by sunpy, and
    therefore only results from 1996 onwards are returned by this client.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.soon)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SRSClient:
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999       SOON ...   SWPC     NOAA
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999       SOON ...   SWPC     NOAA
    <BLANKLINE>
    <BLANKLINE>

    """
    pattern = ('https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_region_summaries/'
               '{{year:4d}}/{{month:2d}}/{{year:4d}}{{month:2d}}{{day:2d}}SRS.txt')

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [("SOON", "Solar Region Summary."),
                                    ("SRS-Table", "Solar Region Summary.")],
                 attrs.Physobs: [('SRS', 'Solar Region Summary.')],
                 attrs.Source: [('SWPC', 'The Space Weather Prediction Center.')],
                 attrs.Provider: [('NOAA', 'The National Oceanic and Atmospheric Administration.')]}
        return adict
