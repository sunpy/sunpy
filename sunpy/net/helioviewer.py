"""
This module provides a wrapper around the Helioviewer API.
"""
from __future__ import absolute_import

#pylint: disable=E1101,F0401,W0231
__author__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
import json
import urllib
import urllib2
import sunpy
from sunpy.time import parse_time

# Helioviewer API URL
__BASE_API_URL__ = "http://helioviewer.org/api/"

def get_data_sources(**kwargs):
    """Returns a structured list of datasources available at Helioviewer.org"""
    params = {"action": "getDataSources"}
    params.update(kwargs)
    
    return _get_json(params)    

def get_closest_image(date, **kwargs):
    """Finds the closest image available for the specified source and date.
    
    For more information on what types of requests are available and the
    expected usage for the response, consult the Helioviewer API documenation:
        http://helioviewer.org/api
    
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
    >>> from sunpy.net import helioviewer
    >>> metadata = helioviewer.get_closest_image('2012/01/01', sourceId=11)
    >>> print(metadata['date'])
    """
    params = {
        "action": "getClosestImage",
        "date": _format_date(date)
    }
    params.update(kwargs)
    
    return _get_json(params)

def get_jp2_image(date, directory=None, **kwargs):
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
        Directory to download JPEG 2000 image to.
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
    out : Returns a map representation of the requested image or a URI if
    "jpip" parameter is set to True.
    
    Examples
    --------
    >>> from sunpy.net import helioviewer
    >>> aia = helioviewer.get_jp2_image('2012/07/03 14:30:00', observatory='SDO', instrument='AIA', detector='AIA', measurement='171')
    >>> aia.show()
    >>>
    >>> data_sources = helioviewer.get_data_sources()
    >>> lasco = helioviewer.get_jp2_image('2012/07/03 14:30:00', sourceId=data_sources['SOHO']['LASCO']['C2']['white-light']['sourceId'])
    """
    params = {
        "action": "getJP2Image",
        "date": _format_date(date)
    }
    params.update(kwargs)
    
    # JPIP URL response
    if 'jpip' in kwargs:
        return _get_json(params)

    filepath = _get_file(params, directory)

    return sunpy.make_map(filepath)

def _get_json(params):
    """Returns a JSON result as a string"""
    response = _request(params).read()
    return json.loads(response)

def _get_file(params, directory=None):
    """Downloads a file and return the filepath to that file"""
    # Query Helioviewer.org
    response = _request(params)
    
    # JPEG 2000 image response
    if directory is None:
        directory = sunpy.config.get('downloads', 'download_dir')
    
    filename = response.info()['Content-Disposition'][22:-1]
    filepath = os.path.join(directory, filename)
    
    f = open(filepath, 'wb')
    f.write(response.read())
    f.close()
    
    return filepath

def _request(params):
    """Sends an API request and returns the result
    
    Parameters
    ----------
    params : dict
        Parameters to send
        
    Returns
    -------
    out : result of request
    """
    response = urllib2.urlopen(__BASE_API_URL__, urllib.urlencode(params))
        
    return response

def _format_date(date):
    """Formats a date for Helioviewer API requests"""
    return parse_time(date).strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + "Z"

if __name__ == "__main__":
    get_closest_image('2012-07-11', 'sdo', 'aia', 'aia', '171')