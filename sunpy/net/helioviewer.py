"""
This module provides a wrapper around the Helioviewer API.
"""
import json
import codecs
import urllib
from pathlib import Path
from collections import OrderedDict

from astropy.utils.decorators import lazyproperty

import sunpy
import sunpy.util.parfive_helpers as parfive
from sunpy.time import parse_time
from sunpy.util.util import partial_key_match
from sunpy.util.xml import xml_to_dict

__all__ = ['HelioviewerClient']


class HelioviewerClient:
    """Helioviewer.org Client"""

    def __init__(self, url="https://api.helioviewer.org/"):
        """
        Parameters
        ----------
        url : `str`
            Default URL points to the Helioviewer API.
        """
        self._api = url

    @lazyproperty
    def data_sources(self):
        """
        We trawl through the return from ``getDataSources`` to create a clean
        dictionary for all available sourceIDs.

        Here is a list of all of them: https://api.helioviewer.org/docs/v2/#appendix_datasources
        """
        data_sources_dict = dict()
        datasources = self.get_data_sources()
        for name, observ in datasources.items():
            # TRACE only has measurements and is thus nested once
            if name == "TRACE":
                for instr, params in observ.items():
                    data_sources_dict[(name, None, None, instr)] = params['sourceId']
            else:
                for inst, detect in observ.items():
                    for wavelength, params in detect.items():
                        if 'sourceId' in params:
                            data_sources_dict[(name, inst, None, wavelength)] = params['sourceId']
                        else:
                            for wave, adict in params.items():
                                data_sources_dict[(name, inst, wavelength, wave)
                                                  ] = adict['sourceId']

        # Sort the output for printing purposes
        return OrderedDict(sorted(data_sources_dict.items(), key=lambda x: x[1]))

    def get_data_sources(self):
        """
        Return a hierarchical dictionary of the available datasources on helioviewer.org.

        This uses ``getDataSources`` from the Helioviewer API.

        Returns
        -------
        out : `dict`
            A dictionary containing meta-information for each data source that Helioviewer supports.
        """
        params = {"action": "getDataSources"}
        return self._get_json(params)

    def get_closest_image(self, date, observatory=None, instrument=None,
                          detector=None, measurement=None, source_id=None):
        """
        Finds the closest image available for the specified source and date.
        **This does not download any file.**

        This uses `getClosestImage <https://api.helioviewer.org/docs/v2/#OfficialClients>`_
        from the Helioviewer API.

        .. note::
            We can use ``observatory`` and ``measurement`` or ``instrument`` and ``measurement`` to
            get the value for source ID which can then be used to get required information.

        Parameters
        ----------
        date : `astropy.time.Time`, `str`
            A `~sunpy.time.parse_time` parsable string or `~astropy.time.Time`
            object for the desired date of the image
        observatory : `str`
            Observatory name
        instrument : `str`
            Instrument name
        detector : `str`
            Detector name
        measurement : `str`
            Measurement name
        source_id : `int`
            ID number for the required instrument/measurement.
            This can be used directly instead of using the previous parameters.

        Returns
        -------
        out : `dict`
            A dictionary containing meta-information for the closest image matched

        Examples
        --------
        >>> from sunpy.net import helioviewer
        >>> client = helioviewer.HelioviewerClient()  # doctest: +REMOTE_DATA
        >>> metadata = client.get_closest_image('2012/01/01', source_id=11)  # doctest: +REMOTE_DATA
        >>> print(metadata['date'])  # doctest: +REMOTE_DATA
        2012-01-01T00:00:07.000
        """
        if source_id is None:
            source_id = self._get_source_id((observatory, instrument, detector, measurement))

        params = {
            "action": "getClosestImage",
            "date": self._format_date(date),
            "sourceId": source_id
        }
        response = self._get_json(params)

        # Cast date string to Time
        response['date'] = parse_time(response['date'])

        return response

    def download_jp2(self, date, progress=True, observatory=None, instrument=None, detector=None,
                     measurement=None, source_id=None, directory=None, overwrite=False):
        """
        Downloads the JPEG 2000 that most closely matches the specified time and
        data source.

        This uses `getJP2Image <https://api.helioviewer.org/docs/v2/#JPEG2000>`_
        from the Helioviewer API.

        .. note::
            We can use ``observatory`` and ``measurement`` or ``instrument`` and ``measurement``
            to get the value for source ID which can then be used to get required information.

        Parameters
        ----------
        date : `astropy.time.Time`, `str`
            A string or `~astropy.time.Time` object for the desired date of the image
        progress : `bool`
            Defaults to True.
            If set to False, disables progress bars seen on terminal when
            downloading files.
        observatory : `str`
            Observatory name
        instrument : `str`
            Instrument name
        measurement : `str`
            Measurement name
        detector : `str`
            Detector name
        source_id : `int`
            ID number for the required instrument/measurement.
            This can be used directly instead of using the previous parameters.
        directory : `str`
            Directory to download JPEG 2000 image to.
        overwrite : `bool`
            Defaults to False.
            If set to True, will overwrite any files with the same name.

        Returns
        -------
        out : `str`
            Returns a filepath to the downloaded JPEG 2000 image.

        Examples
        --------
        >>> import sunpy.map
        >>> from sunpy.net import helioviewer
        >>> hv = helioviewer.HelioviewerClient()  # doctest: +REMOTE_DATA
        >>> filepath = hv.download_jp2('2012/07/03 14:30:00', observatory='SDO',
        ...                            instrument='HMI', detector=None, measurement='continuum')  # doctest: +REMOTE_DATA
        >>> filepath = hv.download_jp2('2012/07/03 14:30:00', observatory='SDO', measurement='continuum')  # doctest: +REMOTE_DATA
        >>> filepath = hv.download_jp2('2012/07/03 14:30:00', instrument='HMI', measurement='continuum')  # doctest: +REMOTE_DATA
        >>> aia = sunpy.map.Map(filepath)  # doctest: +REMOTE_DATA
        >>> aia.peek()  # doctest: +SKIP
        """
        if source_id is None:
            source_id = self._get_source_id((observatory, instrument, detector, measurement))

        params = {
            "action": "getJP2Image",
            "date": self._format_date(date),
            "sourceId": source_id,
        }

        return self._get_file(params, progress=progress, directory=directory, overwrite=overwrite)

    def get_jp2_header(self, date, observatory=None, instrument=None, detector=None, measurement=None, jp2_id=None):
        """
        Get the XML header embedded in a JPEG2000 image. Includes the FITS header as well as a section
        of Helioviewer-specific metadata.

        This uses `getJP2Header <https://api.helioviewer.org/docs/v2/#JPEG2000>`_
        from the Helioviewer API.

        .. note::
            We can use ``observatory`` and ``measurement`` or ``instrument`` and ``measurement``
            to get the value for source ID which can then be used to get required information.

        Parameters
        ----------
        date : `astropy.time.Time`, `str`
            A `~sunpy.time.parse_time` parsable string or `~astropy.time.Time`
            object for the desired date of the image
        observatory : `str`
            Observatory name
        instrument : `str`
            Instrument name
        measurement : `str`
            Measurement name
        detector : `str`
            Detector name
        jp2_id : `int`
            Unique JP2 image identifier.
            This can be used directly instead of using the previous parameters.

        Returns
        -------
        out : `dict`
            Returns a dictionary containing the header information of JPEG 2000 image.
            The returned dictionary may have either one or two keys: *fits* and *helioviewer*.

        Examples
        --------
        >>> from sunpy.net import helioviewer
        >>> hv = helioviewer.HelioviewerClient()  # doctest: +REMOTE_DATA
        >>> header = hv.get_jp2_header('2012/07/03', observatory='SDO',
        ...                            instrument='HMI', detector=None, measurement='continuum')  # doctest: +REMOTE_DATA
        >>> # The key 'fits' can be used to get the fits header information
        >>> fits_header = header['fits']  # doctest: +REMOTE_DATA
        >>> # The keys 'helioviewer' can be used to extract the helioviewer specific metadata.
        >>> helioviewer_meta_data = header['helioviewer']  # doctest: +REMOTE_DATA
        """
        if jp2_id is None:
            jp2_id = self.get_closest_image(
                date, observatory, instrument, detector, measurement)['id']

        params = {
            "action": "getJP2Header",
            "id": jp2_id,
        }

        responses = self._request(params)
        # Reads the output from HTTPResponse object and decodes it.
        responses = responses.read().decode('utf-8')
        return xml_to_dict(responses)['meta']

    def download_png(self, date, image_scale, layers, progress=True,
                     directory=None, overwrite=False, watermark=False,
                     events="", event_labels=False,
                     scale=False, scale_type="earth", scale_x=0, scale_y=0,
                     width=4096, height=4096, x0=0, y0=0,
                     x1=None, y1=None, x2=None, y2=None):
        """
        Downloads the PNG that most closely matches the specified time and
        data source.

        This function is different to `~sunpy.net.helioviewer.HelioviewerClient.download_jp2`.
        Here you get PNG images and return more complex images.

        For example you can return an image that has multiple layers composited together
        from different sources.
        Also mark solar features/events with an associated text label.
        The image can also be cropped to a smaller field of view.

        These parameters are not pre-validated before they are passed to Helioviewer API.
        See https://api.helioviewer.org/docs/v2/#appendix_coordinates for more information about
        what coordinates values you can pass into this function.

        This uses `takeScreenshot <https://api.helioviewer.org/docs/v2/#Screenshots>`_ from the Helioviewer API.

        .. note::

            Parameters ``x1``, ``y1``, ``x2`` and ``y2`` are set to `None`.
            If all 4 are set to values, then keywords: ``width``, ``height``, ``x0``, ``y0`` will be ignored.

        Parameters
        ----------
        date : `astropy.time.Time`, `str`
            A `~sunpy.time.parse_time` parsable string or `~astropy.time.Time` object
            for the desired date of the image
        image_scale : `float`
            The zoom scale of the image in arcseconds per pixel.
            For example, the scale of an AIA image is 0.6.
        layers : `str`
            Image datasource layer/layers to include in the screeshot.
            Each layer string is comma-separated with either:
            "[sourceId,visible,opacity]" or "[obs,inst,det,meas,visible,opacity]".
            Multiple layers are: "[layer1],[layer2],[layer3]".
        progress : `bool`, optional
            Defaults to True.
            If set to False, disables progress bars seen on terminal when
            downloading files.
        events : `str`, optional
            Defaults to an  empty string to indicate no feature/event annotations.
            List feature/event types and FRMs to use to annoate the image.
            Example could be "[AR,HMI_HARP;SPoCA,1]" or "[CH,all,1]"
        event_labels : `bool`, optional
            Defaults to False.
            Annotate each event marker with a text label.
        watermark : `bool`, optional
            Defaults to False.
            Overlay a watermark consisting of a Helioviewer logo and
            the datasource abbreviation(s) and timestamp(s) in the screenshot.
        directory : `str`, optional
            Directory to download JPEG 2000 image to.
        overwrite : bool, optional
            Defaults to False.
            If set to True, will overwrite any files with the same name.
        scale : `bool`, optional
            Defaults to False.
            Overlay an image scale indicator.
        scale_type : `str`, optional
            Defaults to Earth.
            What is the image scale indicator will be.
        scale_x : `int`, optional
            Defaults to 0 (i.e, in the middle)
            Horizontal offset of the image scale indicator in arcseconds with respect
            to the center of the Sun.
        scale_y : `int`, optional
            Defaults to 0 (i.e, in the middle)
            Vertical offset of the image scale indicator in arcseconds with respect
            to the center of the Sun.
        x0 : `float`, optional
            The horizontal offset from the center of the Sun.
        y0 : `float`, optional
            The vertical offset from the center of the Sun.
        width : `int`, optional
            Defaults to 4096.
            Width of the image in pixels.
        height : `int`, optional
            Defaults to 4096.
            Height of the image in pixels.
        x1 : `float`, optional
            Defaults to None
            The offset of the image's left boundary from the center
            of the sun, in arcseconds.
        y1 : `float`, optional
            Defaults to None
            The offset of the image's top boundary from the center
            of the sun, in arcseconds.
        x2 : `float`, optional
            Defaults to None
            The offset of the image's right boundary from the
            center of the sun, in arcseconds.
        y2 : `float`, optional
            Defaults to None
            The offset of the image's bottom boundary from the
            center of the sun, in arcseconds.

        Returns
        -------
        out : `str`
            Returns a filepath to the downloaded PNG image.

        Examples
        --------
        >>> from sunpy.net.helioviewer import HelioviewerClient
        >>> hv = HelioviewerClient()  # doctest: +REMOTE_DATA
        >>> file = hv.download_png('2012/07/16 10:08:00', 2.4,
        ...                        "[SDO,AIA,AIA,171,1,100]",
        ...                        x0=0, y0=0, width=1024, height=1024)   # doctest: +REMOTE_DATA
        >>> file = hv.download_png('2012/07/16 10:08:00', 4.8,
        ...                        "[SDO,AIA,AIA,171,1,100],[SOHO,LASCO,C2,white-light,1,100]",
        ...                        x1=-2800, x2=2800, y1=-2800, y2=2800)   # doctest: +REMOTE_DATA
        """
        params = {
            "action": "takeScreenshot",
            "date": self._format_date(date),
            "imageScale": image_scale,
            "layers": layers,
            "eventLabels": event_labels,
            "events": events,
            "watermark": watermark,
            "scale": scale,
            "scaleType": scale_type,
            "scaleX": scale_x,
            "scaleY": scale_y,
            # Returns the image which we do not want a user to change.
            "display": True
        }

        # We want to enforce that all values of x1, x2, y1, y2 are not None.
        # You can not use both scaling parameters so we try to exclude that here.
        if any(i is None for i in [x1, x2, y1, y2]):
            adict = {"x0": x0, "y0": y0,
                     "width": width, "height": height}
        else:
            adict = {"x1": x1, "x2": x2,
                     "y1": y1, "y2": y2}
        params.update(adict)

        return self._get_file(params, progress=progress, directory=directory, overwrite=overwrite)

    def is_online(self):
        """Returns True if Helioviewer is online and available."""
        try:
            self.get_data_sources()
        except urllib.error.URLError:
            return False

        return True

    def _get_json(self, params):
        """Returns a JSON result as a string."""
        reader = codecs.getreader("utf-8")
        response = self._request(params)
        return json.load(reader(response))

    def _get_file(self, params, progress=True, directory=None, overwrite=False):
        """Downloads a file and return the filepath to that file."""
        if directory is None:
            directory = Path(sunpy.config.get('downloads', 'download_dir'))
        else:
            directory = Path(directory).expanduser().absolute()

        downloader = parfive.Downloader(progress=progress, overwrite=overwrite)

        url = urllib.parse.urljoin(self._api,
                                   "?" + urllib.parse.urlencode(params))
        downloader.enqueue_file(url, path=directory)

        res = downloader.download()

        if len(res) == 1:
            return res[0]
        else:
            return res

    def _request(self, params):
        """
        Sends an API request and returns the result.

        Parameters
        ----------
        params : `dict`
            Parameters to send

        Returns
        -------
        out : result of the request
        """
        response = urllib.request.urlopen(
            self._api, urllib.parse.urlencode(params).encode('utf-8'))
        return response

    def _format_date(self, date):
        """Formats a date for Helioviewer API requests"""
        return parse_time(date).isot + "Z"

    def _get_source_id(self, key):
        """
        Returns source_id based on the key.
        """
        source_id_list = list(partial_key_match(key, self.data_sources))
        if len(source_id_list) != 1:
            raise KeyError(f"The values used: {key} do not correspond to one source_id "
                           f"but {len(source_id_list)} source_id(s)."
                           " Please check the list using HelioviewerClient.data_sources.")
        return source_id_list[0]
