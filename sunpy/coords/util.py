from __future__ import division

__all__ = ['diff_rot']
import numpy as np
import datetime
from sunpy.time import parse_time
from sunpy.wcs import convert_hcc_hg

__author__ = ["Jose Ivan Campos Rozo", "Stuart Mumford"]
__all__ = ['diff_rot']


def diff_rot(ddays, latitude, rot_type='howard', frame_time='sidereal'):
    """
    This function computes the change in longitude over days in degrees.

    Parameters
    -----------
    ddays: float or timedelta
        Number of days to rotate over, or timedelta object.

    latitude: float or array-like
        heliographic coordinate latitude in Degrees.

    rot_type: {'howard' | 'snodgrass' | 'allen'}
        howard: Use values for small magnetic features from Howard et al.
        snodgrass: Use Values from Snodgrass et. al
        allen: Use values from Allen, Astrophysical Quantities, and simplier
                equation.

    frame_time: {'sidereal' | 'synodic'}
        Choose 'type of day' time reference frame.

    Returns:
    -------
    longditude_delta: ndarray
        The change in longitude over days (units=degrees)

    See Also
    --------
    IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro

    Howard rotation:
        http://adsabs.harvard.edu/abs/1990SoPh..130..295H

    A review of rotation parameters (including Snodgrass values):
        http://link.springer.com/article/10.1023%2FA%3A1005226402796

    Examples
    --------
    Default rotation calculation over two days at 30 degrees latitude:
        rotation = diff_rot(2, 30)
    Default rotation over two days for a number of latitudes:
        rotation = diff_rot(2, np.linspace(-70, 70, 20))
    With rotation type 'allen':
        rotation = diff_rot(2, np.linspace(-70, 70, 20), 'allen')
    """

    if not isinstance(ddays, datetime.timedelta):
        delta = datetime.timedelta(days=ddays)

    delta_seconds = delta.total_seconds()
    delta_days = delta_seconds / 24 / 3600

    sin2l = (np.sin(np.deg2rad(latitude)))**2
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

    return np.round(rotation_deg, 4)


def rot_xy(x, y, date=datetime.datetime.utcnow(), **kwargs):
    """ Get a solar rotated position for a given time interval.

    Parameters
    -----------
    x: float or numpy ndarray
        x-locations.

    y: float or numpy ndarray
        y-locations.

    interval: Time interval in seconds; positive (negative) values leads to '
              forward (backward) rotation.

    date: date/time at which the sun position is calculated; can be in any
          format accepted by parse_time. If missing, current date/time is
          assumed.
;
; OUTPUTS:
;       RESULT - A (Mx2) array representing rotated positions in arcsecs,
;                RESULT(*,0) being X position and RESULT(*,1) Y position;
;                where M is number of elements in XX (and in YY).
;                If an error occurs, [-9999,-9999] will be returned.
;
;                If OFFLIMB = 1 after the call, there must be some
;                points rotated to the back of the sun. The points remain
;                visible can be determined by RESULT(INDEX,*), and
;                off-limb points will have the value of (-9999, -9999).
;
; OPTIONAL OUTPUTS:
;       None.
;
; KEYWORDS:
;       DATE    - Date/time at which the sun position is calculated; can
;                 be in any UTC format. If missing, current date/time is
;                 assumed.
;       TSTART  - Date/time to which XX and YY are referred; can be in
;                 any acceptable time format. Must be supplied if
;                 INTERVAL is not passed
;       TEND    - Date/time at which XX and YY will be rotated to; can be
;                 in any acceptable time format. If needed but missing,
;                 current time is assumed

;       OFFLIMB - A named variable indicating whether any rotated
;                 point is off the limb (1) or not (0). When OFFLIMB
;                 is 1, the points still remaining visible (inside the limb)
;                 will be those whose indices are INDEX (below)
;       INDEX   - Indices of XX/YY which remain inside the limb after
;                 rotation. When OFFLIMB becomes 1, the number of
;                 INDEX will be smaller than that of XX or YY. If no
;                 point remains inside the limb, INDEX is set to -1
;       BACK_INDEX  - Indices of XX/YY which were on the disk at TSTART, 
;                 but are no longer visible (i.e. are behind the visible 
;                 solar disk) after rotation. 
;       KEEP    - keep same epoch DATE when rotating; use same P,B0,R values 
;                 both for DATE and DATE+INTERVAL
;       VSTART, VEND = {b0,l0,rsun} = input structure with b0,l0,rsun
;                at start/end of rotation


    Returns:
    -------
    longditude_delta: ndarray
        The change in longitude over days (units=degrees)

    See Also
    --------
    IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/rot_xy.pro

"""
    # must have pairs of co-ordinates
    if np.array(x).shape != np.array(y).shape:
        raise ValueError('Input co-ordinates must have the same shape.')

    # Make sure we have enough time information to perform a solar differential
    # rotation
    if kwargs.get("tend") is not None:
        tend = parse_time(kwargs.get("tend"))
    else:
        tend = datetime.datetime.utcnow()
    if (kwargs.get("tstart") is not None):
        delta = tend - parse_time(kwargs.get("tstart"))

    if kwargs.get("interval") is not None:
        if isinstance(kwargs.get("interval"), datetime.timedelta):
            delta = kwargs.get("interval")
        else:
            delta = datetime.timedelta(seconds=kwargs.get("interval"))
    if ((kwargs.get("tstart") is None) and (kwargs.get("tend") is None)) or \
    ((kwargs.get("tstart") is None) and kwargs.get("interval") is None):
        raise ValueError("You need to specify 'tstart' & 'tend', or " + \
                          "'tstart' and 'interval'")

    # differentially rotate if interval is non-zero
    if delta.total_seconds() > 0.0:
        hgln, hglt = convert_hcc_hg(rsun, b0, l0, x, y)

    return None


def pb0r(date):
    """To calculate the solar P, B0 angles and the semi-diameter.
    http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/pb0r.pro
    """
    pass
