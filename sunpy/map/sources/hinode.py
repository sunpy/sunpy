"""Hinode XRT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm
from sunpy.time import parse_time
from datetime import datetime

class XRTMap(Map):
    """XRT Image Map definition
    
    References
    ----------
    For a description of XRT headers
    """
    #TODO: get a link for the XRT FITS headers
    # Add in some information about the the possible filter wheel measurements
    Map.filter_wheel1_measurements = ["Al_med", "Al_poly", "Be_med",
                                      "Be_thin", "C_poly", "Open"]
    Map.filter_wheel2_measurements = ["Open", "Al_mesh", "Al_thick",
                                      "Be_thick", "Gband", "Ti_poly"]
    @classmethod
    def get_properties(cls, header):
        """Parses XRT image header"""
        properties = Map.get_properties(header)
        # XRT uses DATE_OBS, not date-obs.
        properties["date"] = parse_time(header.get('date_obs', None))

        #TODO: proper exception handling here - report to the user that there is
        # an unexpected value
        fw1 = header.get('EC_FW1_')
        if not(fw1.lower() in [x.lower() for x in cls.filter_wheel1_measurements]):
            pass
        fw2 = header.get('EC_FW2_')
        if not(fw2.lower() in [x.lower() for x in cls.filter_wheel2_measurements]):
            pass

        # All images get the same color table - IDL Red temperature (loadct, 3)
        properties.update({
            "detector": "XRT",
            "instrument": "XRT",
            "observatory": "Hinode",
            "name": "XRT %s-%s " % (fw1.replace('_', ' '),
                                       fw2.replace('_', ' ')),
            "nickname": "XRT",
            "cmap": cm.get_cmap(name='hinodexrt')
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'