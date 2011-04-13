"""Methods relating to the handling of solar data formats

Authors: `Keith Hughitt <keith.hughitt@nasa.gov>`_
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from datetime import datetime
import matplotlib.cm as cm
import matplotlib.colors as colors

#
# Note 2011/04/13: Would it make more sense to instead create a generic class
# which both Map and MapCube inherit from that handles header mapping, etc?
#
def parse_header(header):
    """Parses a FITS, etc image header
    
    Attempts to detect the type of data (e.g. AIA) based on values in the 
    file header and returns a mapped dictionary of of some important values.

    Parameters
    ----------
    header : dict
        A dictionary container the header keywords from the file being read in
    """
    # AIA
    if header['telescop'] == 'SDO/AIA':
        datatype = "aia"
        
    # HMI
    elif header['telescop'] == 'SDO/HMI':
        datatype = "hmi"
        
    # EUVI
    elif header['detector'] == 'EUVI':
        datetype = "euvi"
        
    # COR1/2
    elif (header['detector'] == 'COR1') or (header['detector'] == 'COR2'):
        datatype = "cor"

    # EIT        
    elif header['instrume'] == 'EIT':
        datatype = "eit"
        
    # LASCO
    elif header['instrume'] == 'LASCO':
        datatype = "lasco"
        
    # MDI
    elif header['instrume'] == 'MDI':
        datatype = "mdi"
        
    return _get_norm_header_tags(header, datatype) 

def _get_norm_header_tags(header, type_):
    """Returns a normalized dictionary of header values
    
    A normalized mapping of important header values is created and returned.
    Not all of the header values are used, but instead only those that are
    required for the Map class to function are included. Note that some values
    may be cast to new types in the process.
    
    Parameters
    ----------
    header : dict
        A dictionary container the header keywords from the file being read in
    type_ : str
        A short string describing the type of data being mapped
        
    TODO
    ----
    Rename method to reflect the fact that is now does more than just map
    header values (i.e. cmap, norm)
    
    Returns
    -------
    out : dict
        A new mapped dictionary of useful header values
    """
    date_fmt1 = "%Y-%m-%dT%H:%M:%S.%f"
    date_fmt2 = "%Y-%m-%dT%H:%M:%S.%fZ"
    
    if type_ == "aia":
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "AIA",
            "inst": "AIA",
            "meas": header['wavelnth'],
            "obs": "SDO",
            "name": "AIA %s" % header['wavelnth'],
            "r_sun": header['rsun_obs']
        }
    elif type_ == "hmi":
        return {
            "cmap": cm.gray,
            "norm": None,
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "HMI",
            "inst": "HMI",
            "meas": header['content'].lower(),
            "obs": "SDO",
            "name": "HMI %s" % header['content'].lower(),
            "r_sun": header['rsun_obs']
        }
    elif type_ == "euvi":
        return {
            "cmap": cm.gray,
            "norm": None,
            "date": datetime.strptime(header['date_obs'], date_fmt1),
            "det": "EUVI",
            "inst": "SECCHI",
            "meas": header['wavelnth'],
            "obs": header['obsrvtry'],
            "name": "EUVI %s" % header['wavelnth'],
            "r_sun": header['rsun']
        }
    elif type_ == "cor":
        return {
            "cmap": cm.gray,
            "norm": None,
            "date": datetime.strptime(header['date_obs'], date_fmt1),
            "det": header['detector'],
            "inst": "SECCHI",
            "meas": header['wavelnth'],
            "obs": header['obsrvtry'],
            "name": "SECCHI %s" % header['detector'],
            "r_sun": header['rsun']
        }
    elif type_ == "eit":
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "EIT",
            "inst": "EIT",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "EIT %s" % header['wavelnth'],
            "r_sun": header['solar_r']
        }
    elif type_ == "lasco":
        datestr = "%sT%s" % (header['date_obs'], header['time_obs'])
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(datestr, date_fmt1),
            "det": header['detector'],
            "inst": "LASCO",
            "meas": header['wavelnth'],
            "obs": "SOHO",
            "name": "LASCO %s" % header['detector'],
            "r_sun": None
        }
    elif type_ == "mdi":
        datestr = header['date_obs']
            
        # MDI sometimes has an "60" in seconds field
        if datestr[17:19] == "60":
            datestr = datestr[:17] + "30" + datestr[19:]
        
        # Measurement
        dpcobsr = header['dpc_obsr']
        meas = "Magnetogram" if dpcobsr.find('Mag') != -1 else "Continuum"
        
        return {
            "cmap": None,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(datestr, date_fmt2),
            "det": "MDI",
            "inst": "MDI",
            "meas": meas,
            "obs": "SOHO",
            "name": "MDI %s" % meas,
            "r_sun": header['r_sun']
        }
