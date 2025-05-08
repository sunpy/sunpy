import numpy as np

import astropy.units as u

from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ["AIASynopsisClient"]


class AIASynopsisClient(GenericClient):
    """
    A client for retrieving AIA "synoptic" data from the JSOC.

    This dataset is not synoptic like HMI and MDI Synoptic images which are images of the solar surface reconstructed from many observations over a solar rotation but rather a synopsis of AIA data.

    The AIA synoptic data are calibrated Level 1.5 images with reduced 1k x 1k resolution at regular 2-minute cadence.

    Note
    ----

    This client does not support multiple wavelengths in a single query. If you want to download multiple wavelengths, you need to do it in separate queries.

    References
    ----------
    * `Readme <https://jsoc1.stanford.edu/data/aia/synoptic/README.html>`__

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/1 00:01:00"),
    ...                       a.Instrument.aia, a.Level("1.5s"))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    9 Results from the AIASynopsisClient:
    Source: https://jsoc1.stanford.edu/data/aia/synoptic/
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs  Source Provider Level Wavelength
    ----------------------- ----------------------- ---------- --------- ------ -------- ----- ----------
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S         94
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        131
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        171
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        193
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        211
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        304
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        335
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S       1600
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S       4500
    <BLANKLINE>
    <BLANKLINE>
    """
    pattern = ('https://jsoc1.stanford.edu/data/aia/synoptic/'
               '{{year:4d}}/{{month:2d}}/{{day:2d}}/H{{hour:2d}}00/AIA{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}_{{Wavelength:04d}}.fits')
    required = {a.Time, a.Instrument, a.Level}

    @property
    def info_url(self):
        return "https://jsoc1.stanford.edu/data/aia/synoptic/"

    def search(self, *args, **kwargs):
        """
        Query this client for a list of results.

        Parameters
        ----------
        \\*args: `tuple`
            `sunpy.net.attrs` objects representing the query.
        \\*\\*kwargs: `dict`
            Any extra keywords to refine the search.

        Returns
        -------
        A `QueryResponse` instance containing the query result.
        """
        # TODO: There has to be a better way than repeating the entire search method.
        _, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        if "Wavelength" in matchdict:
            # The scarper uses string matching, so we have to convert the wavelength to a string
            # and remove the trailing zeros to match the pattern on the server.
            matchdict["Wavelength"] = str(np.round(matchdict["Wavelength"].min.to_value(u.AA, u.equivalencies.spectral()))).replace(".0","")
        scraper = Scraper(format=pattern)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        filesmeta = scraper._extract_files_meta(tr, matcher=matchdict)
        filesmeta = sorted(filesmeta, key=lambda k: k['url'])
        metalist = []
        if filesmeta:
            base_time = self.post_search_hook(filesmeta[0], matchdict)['Start Time']
            for i in filesmeta:
                rowdict = self.post_search_hook(i, matchdict)
                # Have to manually deal with sample rate if is it provided
                if sample_rate := matchdict.get("Sample"):
                    sample_rate = u.Quantity(sample_rate[0], u.s)
                    if (rowdict['Start Time'] - base_time).to(u.second) % sample_rate != 0 * u.s:
                        continue
                metalist.append(rowdict)
        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [("AIA", "Data from the Atmospheric Imaging Assembly instrument.")],
            a.Physobs: [
                ("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")
            ],
            a.Source: [("SDO", "The Solar Dynamics Observatory.")],
            a.Provider: [("JSOC", "Joint Science Operations Center at Stanford.")],
            a.Level: [("1.5s", "Level 1.5 data processed for quicker analysis.")],
        }
        return adict
