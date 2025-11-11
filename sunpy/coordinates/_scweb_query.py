"""
Contains helper functions for generating and sending XML payloads to the SSCWeb API
(https://sscweb.gsfc.nasa.gov/) to retrieve spacecraft location data.
"""

import requests
from lxml import etree
from lxml.builder import E

from astropy.time import Time

from sunpy.time import TimeRange

__all__ = ["_create_xml_request", "_send_requests"]

def _create_xml_request(name, time_range, system):
    """
    Create an XML request payload for SSCWeb to fetch location data.

    Following the `API specification <https://sscweb.gsfc.nasa.gov/WebServices/REST/>`__.

    Parameters
    ----------
    name : str
        Name of the satellite or spacecraft.
    time_range : sunpy.time.TimeRange
        The time range for which to fetch location data.

    Returns
    -------
    str
        XML string formatted for SSCWeb API requests.
    """
    if not (isinstance(time_range, TimeRange) or isinstance(time_range, Time)):
        raise ValueError("time_range must be a sunpy.time.TimeRange or astropy.time.Time")

    if isinstance(time_range,Time):
        time_range = TimeRange(time_range[0],time_range[1])

    ns = "http://sscweb.gsfc.nasa.gov/schema"
    xml_request = E.DataRequest(
        {"xmlns": ns},
        E.TimeInterval(
            E.Start(time_range.start.isot),
            E.End(time_range.end.isot),
        ),
        E.Satellites(
            E.Id(str(name).lower()),
            E.ResolutionFactor("1"),
        ),
        E.OutputOptions(
            E.CoordinateOptions(
                E.CoordinateSystem(system.capitalize()),
                E.Component("Lat"),
            ),
            E.CoordinateOptions(
                E.CoordinateSystem(system.capitalize()),
                E.Component("Lon"),
            ),
        ),
    )

    return etree.tostring(xml_request,encoding="unicode")


def _send_requests(xml):
    """
    Send the XML request to `SSCWeb API <https://sscweb.gsfc.nasa.gov/WS/sscr/2/locations>`__ to fetch location data.

    Parameters
    ----------
    xml : str
        XML request payload.

    Returns:
    -------
    requests.Response
        The HTTP response object from the SSCWeb server.
    """
    url = "https://sscweb.gsfc.nasa.gov/WS/sscr/2/locations"

    headers = {
            'Content-Type': 'application/xml',
            'Accept': 'application/xml',
        }
    with requests.Session() as session:
        session.headers.update(headers)
        response = session.post(url, data=xml)
    return response
