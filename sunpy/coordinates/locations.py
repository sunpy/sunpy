import xml.etree.ElementTree as ET

import numpy as np
import requests

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.coordinates import GeocentricSolarEcliptic
from sunpy.time import TimeRange


def _create_xml_request(name, time_range):
    """
    Create an XML request payload for SSCWeb to fetch location data.

    Parameters:
    ----------
    name : str
        Name of the satellite or spacecraft.
    time_range : sunpy.time.TimeRange
        The time range for which to fetch location data.

    Returns:
    -------
    str
        XML string formatted for SSCWeb API requests.
    """
    if not isinstance(time_range, TimeRange):
        raise ValueError("time_range must be a SunPy TimeRange object.")

    # Extract start and end times in ISO format
    start_time = time_range.start.isot
    end_time = time_range.end.isot

    # Root element
    data_request = ET.Element("DataRequest", xmlns="http://sscweb.gsfc.nasa.gov/schema")

    # Description
    ET.SubElement(data_request, "Description")

    # TimeInterval
    time_interval = ET.SubElement(data_request, "TimeInterval")
    ET.SubElement(time_interval, "Start").text = start_time
    ET.SubElement(time_interval, "End").text = end_time

    # Satellites
    satellites = ET.SubElement(data_request, "Satellites")
    ET.SubElement(satellites, "Id").text = name
    ET.SubElement(satellites, "ResolutionFactor").text = "1"

    # OutputOptions
    output_options = ET.SubElement(data_request, "OutputOptions")

    # CoordinateOptions
    coordinate_components = ["Lat", "Lon"]
    for component in coordinate_components:
        coordinate_option = ET.SubElement(output_options, "CoordinateOptions")
        ET.SubElement(coordinate_option, "CoordinateSystem").text = "Gse"
        ET.SubElement(coordinate_option, "Component").text = component

    return ET.tostring(ET.ElementTree(data_request).getroot(), encoding="unicode")

def _send_requests(xml):
    """
    Send the XML request to SSCWeb API to fetch location data.

    Parameters:
    ----------
    xml : str
        XML request payload.

    Returns:
    -------
    requests.Response
        The HTTP response object from the SSCWeb server.
    """
    session = requests.Session()
    url = "https://sscweb.gsfc.nasa.gov/WS/sscr/2/locations"
    headers = {
        'Content-Type': 'application/xml',
        'Accept': 'application/xml',
    }
    session.headers.update(headers)
    response = session.post(url, data=xml)
    return response

def _decode_xml_response(response):
    """
    Decode the XML response from SSCWeb API into a SkyCoord object.

    Parameters:
    ----------
    response : requests.Response
        The HTTP response object containing XML data.

    Returns:
    -------
    astropy.coordinates.SkyCoord
        A SkyCoord object containing latitude, longitude, and timestamps in
        the Geocentric Solar Ecliptic (GSE) frame.
    """
    namespace = {'ns': 'http://sscweb.gsfc.nasa.gov/schema'}
    root = ET.fromstring(response.text)
    coordinates = root.find('.//ns:Coordinates', namespace)

    # Extract latitude and longitude
    latitude_values = np.array([float(lat.text) for lat in coordinates.findall('.//ns:Latitude', namespace)]) * u.deg
    longitude_values = np.array([float(lon.text) for lon in coordinates.findall('.//ns:Longitude', namespace)]) * u.deg

    # Extract time values
    times = root.findall(".//ns:Time", namespace)
    time_values = Time([time.text for time in times])

    return SkyCoord(
        longitude_values,
        latitude_values,
        frame=GeocentricSolarEcliptic,
        obstime=time_values
    )

def get_locations(name, time_range):
    """
    Retrieve the locations of a satellite or spacecraft from SSCWeb API.

    Parameters:
    ----------
    name : str
        Name of the satellite or spacecraft.
    time_range : sunpy.time.TimeRange
        The time range for which to fetch location data.

    Returns:
    -------
    astropy.coordinates.SkyCoord
        A SkyCoord object containing latitude, longitude, and timestamps in
        the GeocentricSolarEcliptic (GSE) coordinate.
    """
    xml = _create_xml_request(name, time_range)
    response = _send_requests(xml)
    return _decode_xml_response(response)
