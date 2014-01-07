"""
This module provides a wrapper around the Helioviewer API.
"""
from __future__ import absolute_import

#pylint: disable=E1101,F0401,W0231

__author__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
import urllib
import urllib2

import json

import sunpy
from sunpy.time import parse_time
from sunpy.util.net import download_fileobj

__all__ = ['HelioviewerClient']

class HelioviewerClient:
    """Helioviewer.org Client"""
    def __init__(self, url="http://helioviewer.org/api/"):
        self._api = url

    def get_data_sources(self, **kwargs):
        """Returns a structured list of datasources available at Helioviewer.org"""
        params = {"action": "getDataSources"}
        params.update(kwargs)

        return self._get_json(params)

    def get_closest_image(self, date, **kwargs):
        """Finds the closest image available for the specified source and date.

        For more information on what types of requests are available and the
        expected usage for the response, consult the Helioviewer
        API documenation: http://helioviewer.org/api

        Parameters
        ----------
        date : datetime, string
            A string or datetime object for the desired date of the image
        observatory : string
            (Optional) Observatory name
        instrument : string
            (Optional) instrument name
        detector : string
            (Optional) detector name
        measurement : string
            (Optional) measurement name
        sourceId : int
            (Optional) data source id

        Returns
        -------
        out : dict
            A dictionary containing metainformation for the closest image matched

        Examples
        --------
        >>> from sunpy.net import HelioviewerClient

        >>> client = HelioviewerClient()
        >>> metadata = client.get_closest_image('2012/01/01', sourceId=11)
        >>> print(metadata['date'])
        """
        params = {
            "action": "getClosestImage",
            "date": self._format_date(date)
        }
        params.update(kwargs)

        response = self._get_json(params)

        # Cast date string to DateTime
        response['date'] = parse_time(response['date'])

        return response

    def download_jp2(self, date, directory=None, overwrite=False, **kwargs):
        """
        Downloads the JPEG 2000 that most closely matches the specified time and
        data source.

        The data source may be specified either using it's sourceId from the
        get_data_sources query, or a combination of observatory, instrument,
        detector and measurement.

        Parameters
        ----------
        date : datetime, string
            A string or datetime object for the desired date of the image
        directory : string
            (Optional) Directory to download JPEG 2000 image to.
        observatory : string
            (Optional) Observatory name
        instrument : string
            (Optional) instrument name
        detector : string
            (Optional) detector name
        measurement : string
            (Optional) measurement name
        sourceId : int
            (Optional) data source id
        jpip : bool
            (Optional) Returns a JPIP URI if set to True

        Returns
        -------
        out : string
            Returns a filepath to the downloaded JPEG 2000 image or a URL if
            the "jpip" parameter is set to True.

        Examples
        --------
        >>> import sunpy
        >>> from sunpy.net import helioviewer
        >>> hv = helioviewer.HelioviewerClient()
        >>> filepath = hv.download_jp2('2012/07/03 14:30:00', observatory='SDO', instrument='AIA', detector='AIA', measurement='171')
        >>> aia = sunpy.make_map(filepath)
        >>> aia.show()

        >>> data_sources = hv.get_data_sources()
        >>> hv.download_jp2('2012/07/03 14:30:00', sourceId=data_sources['SOHO']['LASCO']['C2']['white-light']['sourceId'])
        """
        params = {
            "action": "getJP2Image",
            "date": self._format_date(date)
        }
        params.update(kwargs)

        # JPIP URL response
        if 'jpip' in kwargs:
            return self._get_json(params)

        return self._get_file(params, directory, overwrite=overwrite)

    def download_png(self, date, image_scale, layers, directory=None,
                     overwrite=False, **kwargs):
        """Downloads a PNG image using data from Helioviewer.org.

        Returns a single image containing all layers/image types requested.
        If an image is not available for the date requested the closest
        available image is returned. The region to be included in the
        image may be specified using either the top-left and bottom-right
        coordinates in arc-seconds, or a center point in arc-seconds and a
        width and height in pixels. See the Helioviewer.org API Coordinates
        Appendix for more infomration about working with coordinates in
        Helioviewer.org.

        Parameters
        ----------
        date : datetime, string
            A string or datetime object for the desired date of the image
        image_scale : float
            The zoom scale of the image. Default scales that can be used are
            0.6, 1.2, 2.4, and so on, increasing or decreasing by a factor
            of 2. The full-res scale of an AIA image is 0.6.
        layers : string
            Each layer string is comma-separated with these values, e.g.:
            "[sourceId,visible,opacity]" or "[obs,inst,det,meas,visible,opacity]"
            Mulitple layer string are by commas: "[layer1],[layer2],[layer3]"
        directory : string
            (Optional)  Directory to download JPEG 2000 image to.
        x1 : float
            (Optional) The offset of the image's left boundary from the center
            of the sun, in arcseconds.
        y1 : float
            (Optional) The offset of the image's top boundary from the center
            of the sun, in arcseconds.
        x2 : float
            (Optional) The offset of the image's right boundary from the
            center of the sun, in arcseconds.
        y2 : float
            (Optional) The offset of the image's bottom boundary from the
            center of the sun, in arcseconds.
        x0 : float
            (Optional) The horizontal offset from the center of the Sun.
        y0 : float
            (Optional) The vertical offset from the center of the Sun.
        width : int
            (Optional) Width of the image in pixels (Maximum: 1920).
        height : int
            (Optional) Height of the image in pixels (Maximum: 1200).
        watermark
            (Optional) Whether or not the include the timestamps and the
            Helioviewer.org logo in the image (Default=True).

        Returns
        -------
        out : string
            filepath to the PNG image

        Examples
        --------
        >>> from sunpy.net.helioviewer import HelioviewerClient
        >>> hv = HelioviewerClient()
        >>> hv.download_png('2012/07/16 10:08:00', 2.4, "[SDO,AIA,AIA,171,1,100]", x0=0, y0=0, width=1024, height=1024)
        '/home/user/sunpy/data/2012_07_16_10_08_00_AIA_171.png
        >>> hv.download_png('2012/07/16 10:08:00', 4.8, "[SDO,AIA,AIA,171,1,100],[SOHO,LASCO,C2,white-light,1,100]", x1=-2800, x2=2800, y1=-2800, y2=2800, directory='~/Desktop')
        '/home/user/Desktop/2012_07_16_10_08_00_AIA_171__LASCO_C2.png'
        """
        params = {
            "action": "takeScreenshot",
            "date": self._format_date(date),
            "imageScale": image_scale,
            "layers": layers,
            "display": True
        }
        params.update(kwargs)

        return self._get_file(params, directory, overwrite=overwrite)

    def is_online(self):
        """Returns True if Helioviewer is online and available"""
        try:
            self.get_data_sources()
        except urllib2.URLError:
            return False

        return True

    def _get_json(self, params):
        """Returns a JSON result as a string"""
        response = self._request(params).read()
        return json.loads(response)

    def _get_file(self, params, directory=None, overwrite=False):
        """Downloads a file and return the filepath to that file"""
        # Query Helioviewer.org
        if directory is None:
            directory = sunpy.config.get('downloads', 'download_dir')
        else:
            directory = os.path.abspath(os.path.expanduser(directory))

        response = self._request(params)
        try:
            filepath = download_fileobj(response, directory, overwrite=overwrite)
        finally:
            response.close()

        return filepath

    def _request(self, params):
        """Sends an API request and returns the result

        Parameters
        ----------
        params : dict
            Parameters to send

        Returns
        -------
        out : result of request
        """
        response = urllib2.urlopen(self._api, urllib.urlencode(params))

        return response

    def _format_date(self, date):
        """Formats a date for Helioviewer API requests"""
        return parse_time(date).strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + "Z"
