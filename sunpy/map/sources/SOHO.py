"""SOHO Map subclass definitions"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.BaseMap import BaseMap
from datetime import datetime

class EITMap(BaseMap):
    """EIT Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        date_format = "%Y-%m-%dT%H:%M:%S.%fZ"
        
        properties = BaseMap.get_properties()
        properties.update({
            "date": datetime.strptime(header['date_obs'], date_format),
            "det": "EIT",
            "inst": "EIT",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "EIT %s" % header['wavelnth'],
            "r_sun": header['solar_r']
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'

class LASCOMap(BaseMap):
    """LASCO Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        datestr = "%sT%s" % (header['date_obs'], header['time_obs'])

        properties = BaseMap.get_properties()
        properties.update({
            "date": datetime.strptime(datestr, "%Y/%m/%dT%H:%M:%S.%f"),
            "det": header['detector'],
            "inst": "LASCO",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "LASCO %s" % header['detector'],
            "r_sun": None
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an LASCO image"""
        return header.get('instrume') == 'LASCO'
        
class MDIMap(BaseMap):
    """MDI Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
        
    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        # MDI sometimes has an "60" in seconds field
        datestr = header['date_obs']

        if datestr[17:19] == "60":
            datestr = datestr[:17] + "30" + datestr[19:]
        
        # Determine measurement
        dpcobsr = header['dpc_obsr']
        meas = "magnetogram" if dpcobsr.find('Mag') != -1 else "continuum"
        
        properties = BaseMap.get_properties()        
        properties.update({
            "date": datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S.%fZ"),
            "det": "MDI",
            "inst": "MDI",
            "meas": meas,
            "obs": "SOHO",
            "name": "MDI %s" % meas,
            "r_sun": header['radius']
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an MDI image"""
        return header.get('instrume') == 'MDI'

