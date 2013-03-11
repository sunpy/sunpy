__author__ = ["Jose Ivan Campos Rozo"]

__all__ = ['diff_rot']
import numpy as np
import datetime

def diff_rot(ddays,latitude,rot_type=None):
    """
    This function computes the change in longitude over days in degrees.

    Parameters
    -----------
    ddays: float or timedelta
        Number of days that I want to rotate.
        
    latitude: float or array-like
        heliographic coordinate latitude in Degrees.
        
    rot_type: {None | 'howard' | 'synodic' | 'sidereal' | 'allen'}
        howard: Use values for small magnetic features from Howard et al.
        synodic: Use synodic rotation rate.
        sidereal: Use sidereal rotation rate.
        allen: Use values from Allen, Astrophysical Quantities

    Returns:
    -------
    longditude_delta: ndarray            
        The change in longitude over days (units=degrees)
    
    See Also
    --------
    IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro
    
    Howard rotation:
        http://labs.adsabs.harvard.edu/adsabs/abs/1990SoPh..130..295H
    
    Examples
    --------
    Default rotation calculation over two days at 30 degrees latitude:
        rotation = diff_rot(2, 30)
    Defult rotation over two days for a number of latitudes:
        rotation = diff_rot(2, np.linspace(-70, 70, 20))
    With rotation type 'allen':
        rotation = diff_rot(2, np.linspace(-70, 70, 20), 'allen')
    """
    if not isinstance(ddays,datetime.timedelta):
        delta = datetime.timedelta(days=ddays)
    
    sin2l = (np.sin(np.deg2rad(latitude)))**2
    sin4l = sin2l**2
    
    if rot_type in [None, 'howard', 'sidereal', 'synodic']:
	    rotation = 1. * 10**-6 * delta.days * (2.894 - 0.428 * sin2l -
                                    0.37 * sin4l) * 24. * 3600. / np.deg2rad(1)
	    if rot_type == 'synodic':
		    rotation = rotation - 0.9856 * delta.days
		 
    elif rot_type == 'allen':
        rotation = delta.days * (14.44 - (3.0 * sin2l))
    
    else:
        raise ValueError("""rot_type must equal one of 
        {None | 'howard' | 'synodic' | 'sidereal' | 'allen'}""")

    return np.round(rotation,4)