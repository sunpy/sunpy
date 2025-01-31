"""
Contains helper functions for generating and sending XML payloads to the SSCWeb API
(https://sscweb.gsfc.nasa.gov/) to retrieve spacecraft location data.
"""

import xml.etree.ElementTree as ET

import requests

from sunpy.time import TimeRange


def _create_xml_request(name, time_range,system):
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
    start_time = time_range.start.isot
    end_time = time_range.end.isot
    data_request = ET.Element("DataRequest", xmlns="http://sscweb.gsfc.nasa.gov/schema")
    ET.SubElement(data_request, "Description")
    time_interval = ET.SubElement(data_request, "TimeInterval")
    ET.SubElement(time_interval, "Start").text = start_time
    ET.SubElement(time_interval, "End").text = end_time
    satellites = ET.SubElement(data_request, "Satellites")
    ET.SubElement(satellites, "Id").text = name
    ET.SubElement(satellites, "ResolutionFactor").text = "1"
    output_options = ET.SubElement(data_request, "OutputOptions")
    coordinate_components = ["Lat", "Lon"]
    for component in coordinate_components:
        coordinate_option = ET.SubElement(output_options, "CoordinateOptions")
        ET.SubElement(coordinate_option, "CoordinateSystem").text = system
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
