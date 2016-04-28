from __future__ import absolute_import
from __future__ import division


import os
import copy
from sunpy.extern.six.moves import urllib
from collections import OrderedDict
import tempfile

import datetime
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Longitude, Latitude

from sunpy.time import parse_time, TimeRange
from sunpy import sun
from sunpy.io.fits import fits

__all__ = ['download_weekly_pointing_file', 'get_detector_sun_angles_for_time',
            'get_detector_sun_angles_for_date','plot_detector_sun_angles',
             'met_to_utc']

def download_weekly_pointing_file(date):
    """
    Downloads the FERMI/LAT weekly pointing file corresponding to the specified date.
    This file contains 1 minute cadence data on the spacecraft pointing, useful for
    calculating detector angles.

    Parameters
    ----------

    date : `datetime.datetime`
        A datetime object or other date format understood by the parse_time function.
    """

    date = parse_time(date)
    # use a temp directory to hold the file
    tmp_dir = tempfile.mkdtemp()
    # use Fermi data server to access weekly LAT pointing file.
    base_url = 'ftp://legacy.gsfc.nasa.gov/FTP/glast/data/lat/weekly/spacecraft/'
    fbasename = 'lat_spacecraft_weekly_w'
    
    # find out which file to get based on date
    # earliest full file in the FERMI server is for mission week 10,
    # beginning 2008 August 7.
    weekly_file_start = parse_time('2008-08-07')
    base_week = 10

    # find out which mission week corresponds to date
    time_diff = date-weekly_file_start
    weekdiff = time_diff.days//7
    week = weekdiff + base_week
   # weekstr = ('%03.0f' % week)
    weekstr = '{:03.0f}'.format(week)

    # construct the full url for the weekly pointing file
    full_fname = fbasename + weekstr + '_p202_v001.fits'
    pointing_file_url = base_url + full_fname

    #try to download the file from the FTP site
    try:
        resp = urllib.request.urlopen(pointing_file_url)
        exists = True
    except:
        urllib.error.HTTPError
        exists = False

    # if no matches at all were found, then the pointing file doesn't exist
    if not exists:
        raise ValueError('No Fermi pointing files found for given date!')

    # download the file
    destination=os.path.join(tmp_dir,full_fname)
    urllib.request.urlretrieve(pointing_file_url,destination)

    # return the location of the downloaded file
    return destination


def get_detector_sun_angles_for_time(time, file):
    """
    get the GBM detector angles vs the sun for a single time.

    Parameters
    ----------

    time : `datetime.datetime`
        A datetime object or other time format understood by the parse_time function.
    file : `str`
        A filepath to a Fermi/LAT weekly pointing file (e.g. as obtained by the
        download_weekly_pointing_file function).
    """

    time = parse_time(time)
    scx, scz, tt = get_scx_scz_at_time(time, file)
    # retrieve the detector angle information in spacecraft coordinates
    detectors = nai_detector_angles()

    # get the detector pointings in RA/DEC given the input spacecraft x and z
    # axes
    detector_radecs = nai_detector_radecs(detectors, scx, scz, tt)

    # this gets the sun position with RA in hours in decimal format (e.g. 4.3).
    # DEC is already in degrees
    sunpos_ra_not_in_deg = [sun.sun.apparent_rightascension(time),
                            sun.sun.apparent_declination(time)]
    # now Sun position with RA in degrees
    sun_pos = [sunpos_ra_not_in_deg[0].to('deg'), sunpos_ra_not_in_deg[1]]
    # sun_pos = [(sunpos_ra_not_in_deg[0] / 24) * 360., sunpos_ra_not_in_deg[1]]
    # now get the angle between each detector and the Sun
    detector_to_sun_angles = (get_detector_separation_angles(detector_radecs, sun_pos))

    return detector_to_sun_angles


def get_detector_sun_angles_for_date(date, file):
    """
    get the GBM detector angles vs the sun as a function of time for a given date

    Parameters
    ----------

    date : `datetime.datetime`
        A datetime object or other date format understood by the parse_time function.
    file : `str`
        A filepath to a Fermi/LAT weekly pointing file (e.g. as obtained by the
        download_weekly_pointing_file function).
    """

    date = parse_time(date)
    tran = TimeRange(date, date + datetime.timedelta(days=1))
    scx, scz, times = get_scx_scz_in_timerange(tran, file)

    # retrieve the detector angle information in spacecraft coordinates
    detectors = nai_detector_angles()

    detector_to_sun_angles = []
    # get the detector vs Sun angles for each t and store in a list of
    # dictionaries.
    for i in range(len(scx)):
        detector_radecs = nai_detector_radecs(detectors, scx[i], scz[i], times[i])

        # this gets the sun position with RA in hours in decimal format
        # (e.g. 4.3). DEC is already in degrees
        sunpos_ra_not_in_deg = [sun.sun.apparent_rightascension(times[i]),
                                sun.sun.apparent_declination(times[i])]
        # now Sun position with RA in degrees
        sun_pos = [sunpos_ra_not_in_deg[0].to('deg'), sunpos_ra_not_in_deg[1]]
        # now get the angle between each detector and the Sun
        detector_to_sun_angles.append(get_detector_separation_angles(detector_radecs, sun_pos))

    # slice the list of dictionaries to get the angles for each detector in a
    # list form
    angles = OrderedDict()
    key_list = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','time']
    for i in range(13):
        if not key_list[i] == 'time':
            angles[key_list[i]] = [item[key_list[i]].value
                                   for item in detector_to_sun_angles] * u.deg
        else:
            angles[key_list[i]] = [item[key_list[i]]
                                   for item in detector_to_sun_angles]

    return angles


def plot_detector_sun_angles(angles):
    """
    Plots the Fermi/GBM detector angles as a function of time.

    Parameters
    ----------

    angles : `dict`
        A dictionary containing the Fermi/GBM detector angle information as a
        function of time. Obtained from the get_detector_separation_angles
        function.
    """

    # make a plot showing the angles vs time
    figure = plt.figure(1)
    for n in angles.keys():
        if not n == 'time':
            plt.plot(angles['time'],angles[n].value,
                     label = '{lab} ({val})' .format(lab=n,
                     val = str(np.mean(angles[n].value))[0:5]))
    plt.ylim(180,0)
    plt.ylabel('angle (degrees)')
    plt.xlabel('Start time: ' + angles['time'][0].isoformat())
    plt.title('Detector pointing angle from Sun')
    plt.legend(fontsize=10)
    figure.autofmt_xdate()
    plt.show()


def get_scx_scz_at_time(time, file):
    """
    Read a downloaded FERMI weekly pointing file and extract scx, scz for a
    single time.

    Parameters
    ----------

    time : `datetime.datetime`
        A datetime object or other time format understood by the parse_time function.
    file : `str`
        A filepath to a Fermi/LAT weekly pointing file (e.g. as obtained by the
         download_weekly_pointing_file function).
    """

    time = parse_time(time)
    hdulist = fits.open(file)
    timesinutc = []
    for tim in hdulist[1].data['START']:
        timesinutc.append(met_to_utc(tim))
    ind = np.searchsorted(timesinutc, time)

    scx_radec = (Longitude(hdulist[1].data['RA_SCX'][ind]*u.deg),
                 Latitude(hdulist[1].data['DEC_SCX'][ind]*u.deg))
    scz_radec = (Longitude(hdulist[1].data['RA_SCZ'][ind]*u.deg),
                 Latitude(hdulist[1].data['DEC_SCZ'][ind]*u.deg))

    return scx_radec, scz_radec, timesinutc[ind]


def get_scx_scz_in_timerange(timerange, file):
    """
    Read a downloaded FERMI weekly pointing file and extract scx, scz for a
    timerange.

    Parameters
    ----------

    timerange : `sunpy.time.TimeRange`
        A SunPy TimeRange object.
    file : `str`
        A filepath to a Fermi/LAT weekly pointing file (e.g. as obtained by the
        download_weekly_pointing_file function).
    """

    hdulist = fits.open(file)
    timesinutc = []
    for tim in hdulist[1].data['START']:
        timesinutc.append(met_to_utc(tim))

    startind = np.searchsorted(timesinutc, timerange.start)
    endind = np.searchsorted(timesinutc, timerange.end)

    scx_radec = []
    scz_radec = []
    for i in range(startind, endind):
        scx_radec.append((Longitude(hdulist[1].data['RA_SCX'][i]*u.deg),
                          Latitude(hdulist[1].data['DEC_SCX'][i]*u.deg)))
        scz_radec.append((Longitude(hdulist[1].data['RA_SCZ'][i]*u.deg),
                          Latitude(hdulist[1].data['DEC_SCZ'][i]*u.deg)))
    return scx_radec, scz_radec, timesinutc[startind:endind]


def nai_detector_angles():
    """
    Returns the dictionary of Fermi/GBM NAI detector zenith and azimuth angles,
    in spacecraft coordinates. zenith angle is measured from +z (along the LAT
    boresight), azimuth is measured from +x. see Meegan et al. (2009) for
    details and detector angles.
    """

    # angles listed as [azimuth, zenith]
    detectors = {'n0': [45.89*u.deg, 20.58*u.deg],
                 'n1': [45.11*u.deg, 45.31*u.deg],
                 'n2': [58.44*u.deg, 90.21*u.deg],
                 'n3': [314.87*u.deg, 45.24*u.deg],
                 'n4': [303.15*u.deg, 90.27*u.deg],
                 'n5': [3.35*u.deg, 89.79*u.deg],
                 'n6': [224.93*u.deg, 20.43*u.deg],
                 'n7': [224.62*u.deg, 46.18*u.deg],
                 'n8': [236.61*u.deg, 89.97*u.deg],
                 'n9': [135.19*u.deg, 45.55*u.deg],
                 'n10': [123.73*u.deg, 90.42*u.deg],
                 'n11': [183.74*u.deg, 90.32*u.deg]
                 }

    return detectors


def nai_detector_radecs(detectors, scx, scz, time):
    """
    calculates the RA/DEC for each NaI detector given spacecraft z and x RA/DEC
    positions.

    NB: This routine is based on code found in GTBURST, originally written by
    Dr Giacamo Vianello for the Fermi Science Tools.

    Parameters
    ----------

    detectors : `dict`
        A dictionary containing the Fermi/GBM detector pointing angles relative
        to the spacecraft axes. Obtained from the nai_detector_angles function.
    scx : array-like
        Two-element tuple containing the RA/DEC information of the Fermi
        spacecraft X-axis
    scz : array-like
        Two-element tuple containing the RA/DEC information of the Fermi
        spacecraft Z-axis
    time : `datetime.datetime`
        The time corresponding to the input scx and scz values, in a format
        understandable by parse_time.

    Returns
    -------
    `dict`
        A dictionary containing the RA/DEC for each Fermi/GBM NaI detector at
        the given input time.
    """

    scx_vector = (np.array([np.cos(scx[0].to('rad').value)*np.cos(scx[1].to('rad').value),
                  np.sin(scx[0].to('rad').value)*np.cos(scx[1].to('rad').value),
                  np.sin(scx[1].to('rad').value)]))

    scz_vector = (np.array([np.cos(scz[0].to('rad').value)*np.cos(scz[1].to('rad').value),
                  np.sin(scz[0].to('rad').value)*np.cos(scz[1].to('rad').value),
                  np.sin(scz[1].to('rad').value)]))

    # For each detector, do the rotation depending on the detector zenith and
    # azimuth angles.
    detector_radecs = copy.deepcopy(detectors)
    for l, d in detectors.items():
        phi = d[0].value
        theta = d[1].value

        # rotate about spacecraft z-axis first
        vx_primed = rotate_vector(scx_vector, scz_vector, np.deg2rad(phi))

        # now find spacecraft y-axis using cross product
        vy_primed = np.cross(scz_vector, vx_primed)

        # do the second part of the rotation around vy
        vz_primed = rotate_vector(scz_vector, vy_primed, np.deg2rad(theta))

        # now we should be pointing at the new RA/DEC.
        ra = Longitude(np.degrees(np.arctan2(vz_primed[1], vz_primed[0])) * u.deg)
        dec = Latitude(np.degrees(np.arcsin(vz_primed[2])) * u.deg)

        # save the RA/DEC in a dictionary
        detector_radecs[l] = [ra, dec]

    detector_radecs['time'] = time
    return detector_radecs


def rotate_vector(vector, axis, theta):
    """
    The Euler-Rodrigues formula for rotating vectors.

    Parameters
    ----------
    vector : `numpy.ndarray`
          a three-element vector to be rotated
    axis : `numpy.ndarray`
          the the-element vector to rotate around
    theta : `float`
          the angle (in radians) by which to rotate vector around axis

    Reference
    ---------
    http://en.wikipedia.org/wiki/Euler-Rodrigues_parameters#Rotation_angle_and_rotation_axis
    """

    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)

    rot_matrix = np.array([[a*a+b*b-c*c-d*d, 2*(b*c+a*d),2*(b*d-a*c)],
                 [2*(b*c-a*d), a*a+c*c-b*b-d*d, 2*(c*d+a*b)],
                 [2*(b*d+a*c), 2*(c*d-a*b), a*a+d*d-b*b-c*c]])

    return np.dot(rot_matrix, vector)


def get_detector_separation_angles(detector_radecs, sunpos):
    """
    Finds the separation angle between the Sun and each NaI detector,
    given a dictionary of detector RA/DECs.

    Parameters
    ----------
    detector_radecs : `dict`
            the RA/DEC for each NaI detector as AstroPy quantities. Obtained
            from the fermi.nai_detector_radecs function
    sunpos : `list`
            Two-element list containing the RA/DEC of the Sun position as
            AstroPy Quantities, e.g. [<Longitude 73.94 deg>,
            <Latitude 22.66 deg>]
    
    """
    angles = copy.deepcopy(detector_radecs)
    for l, d in detector_radecs.items():
        if not l == 'time':
            angle = separation_angle(d, sunpos)
            angles[l] = angle

    return angles


def separation_angle(radec1, radec2):
    """
    Use the law of spherical cosines to calculate the separation angle
    between two RA/DEC positions.

    Parameters
    ----------
    radec1 : `list`
           A two-element list containing an RA/DEC position as AstroPy Quantities,
           e.g. [<Longitude 73.94 deg>, <Latitude 22.66 deg>]
    radec2 : `list`
           A two-element list containing an RA/DEC position as AstroPy Quantities,
           e.g. [<Longitude 73.94 deg>, <Latitude 22.66 deg>]

    """

    cosine_of_angle = (( np.cos( ((90*u.deg) - radec1[1].to('degree')).to('rad')) *
                        np.cos( (90*u.deg - radec2[1].to('degree')).to('rad')) )
                        + ( np.sin( ((90*u.deg) - radec1[1].to('degree')).to('rad') ) *
                        np.sin( ((90*u.deg) - radec2[1].to('degree')).to('rad') ) *
                      np.cos( (radec1[0].to('degree') - radec2[0].to('degree')).to('rad')   ) ) )

    angle = (np.arccos(cosine_of_angle)).to('degree')

    return angle


def met_to_utc(timeinsec):
    """
    Converts Fermi Mission Elapsed Time (MET) in seconds to a datetime object.

    Parameters
    ----------
    timeinsec : `float`
        time in seconds since 00:00 UT on 1st January 2001 (the Fermi MET
        format)

    Returns
    -------
    `datetime.datetime`
        The input Fermi Mission Elapsed Time converted to a datetime object.

    """
    # Times for GBM are in Mission Elapsed Time (MET).
    # The reference time for this is 2001-Jan-01 00:00.
    met_ref_time = parse_time('2001-01-01 00:00')
    offset_from_utc = (met_ref_time - parse_time('1979-01-01 00:00')).total_seconds()
    time_in_utc = parse_time(timeinsec + offset_from_utc)

    return time_in_utc

def utc_to_met(time_ut):
    """
    Converts a UT (in datetime format) to a Fermi Mission Elapsed Time (MET) float.

    Parameters
    ----------
    time_ut : 'datetime.datetime'
        A datetime object in UT

    Returns
    -------
    'float'
        The Fermi Mission Elapsed Time corresponding to the input UT

    """
    met_ref_time = parse_time('2001-01-01 00:00')
    ut_seconds = (time_ut - parse_time('1979-01-01')).total_seconds()
    offset_from_utc = (met_ref_time - parse_time('1979-01-01 00:00')).total_seconds()
    fermi_met = ut_seconds - offset_from_utc

    return fermi_met
