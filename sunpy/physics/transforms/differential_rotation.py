from __future__ import division

__all__ = ['diff_rot']
import numpy as np
import datetime
from astropy import units as u
from astropy.coordinates import Longitude

__author__ = ["Jose Ivan Campos Rozo","Stuart Mumford"]
__all__ = ['diff_rot']


def diff_rot(ddays, latitude, rot_type='howard', frame_time='sidereal'):
    """
    This function computes the change in longitude over days in degrees.

    Parameters
    -----------
    ddays: float or timedelta
        Number of days to rotate over, or timedelta object.

    latitude: '~astropy.units.Quantity` instance
        heliographic coordinate latitude in Degrees.

    rot_type: {'howard' | 'snodgrass' | 'allen'}
        howard: Use values for small magnetic features from Howard et al.
        snodgrass: Use Values from Snodgrass et. al
        allen: Use values from Allen, Astrophysical Quantities, and simplier equation.

    frame_time: {'sidereal' | 'synodic'}
        Choose 'type of day' time reference frame.

    Returns
    -------
    longditude_delta: astropy.units.Quantity     
        The change in longitude over days (units=degrees)

    Notes
    -----
    * IDL code equavalent: http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro
    * Howard rotation: http://adsabs.harvard.edu/abs/1990SoPh..130..295H
    * A review of rotation parameters (including Snodgrass values): http://link.springer.com/article/10.1023%2FA%3A1005226402796

    Examples
    --------
    Default rotation calculation over two days at 30 degrees latitude:
    
    >>> rotation = diff_rot(2, 30 * u.deg)
    
    Default rotation over two days for a number of latitudes:
    
    >>> rotation = diff_rot(2, np.linspace(-70, 70, 20) * u.deg)
    
    With rotation type 'allen':
    
    >>> rotation = diff_rot(2, np.linspace(-70, 70, 20) * u.deg, 'allen')
    """

    if not isinstance(ddays,datetime.timedelta):
        delta = datetime.timedelta(days=ddays)

    if not isinstance(latitude, u.Quantity):
	raise TypeError("Expecting astropy Quantity")

    latitude = latitude.to(u.deg)
    delta_seconds = delta.total_seconds()
    delta_days = delta_seconds / 24 / 3600
    
    sin2l = (np.sin(latitude))**2
    sin4l = sin2l**2

    rot_params = {'howard': [2.894, -0.428, -0.370],
                  'snodgrass': [2.851, -0.343, -0.474]
                  }

    if rot_type not in ['howard', 'allen', 'snodgrass']:
        raise ValueError("""rot_type must equal one of
                        { 'howard' | 'allen' | 'snodgrass' }""")

    elif rot_type == 'allen':
        rotation_deg = delta_days * (14.44 - (3.0 * sin2l))

    else:
        A, B, C = rot_params[rot_type]

        #This is in micro-radians / sec
        rotation_rate = A + B * sin2l + C * sin4l
        rotation_deg = rotation_rate * 1e-6  * delta_seconds / np.deg2rad(1)

    if frame_time == 'synodic':
        rotation_deg -= 0.9856 * delta_days
    
    return Longitude((np.round(rotation_deg,4)),u.deg)
