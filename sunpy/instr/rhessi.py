# -*- coding: utf-8 -*-
"""
    Provides programs to process and analyze RHESSI data.

    .. warning:: This module is in development.

"""
from __future__ import absolute_import, print_function

import csv
from datetime import datetime
from datetime import timedelta

import numpy as np

from astropy.io import fits
from astropy import units as u

import sunpy.map
import sunpy.sun.constants

from sunpy.time import TimeRange, parse_time
from sunpy.sun.sun import solar_semidiameter_angular_size
from sunpy.sun.sun import sunearth_distance

from sunpy.extern.six.moves import urllib

__all__ = ['get_obssumm_dbase_file', 'parse_obssumm_dbase_file',
           'get_obssum_filename', 'get_obssumm_file', 'parse_obssumm_file',
           'backprojection']

# Measured fixed grid parameters
grid_pitch = (4.52467, 7.85160, 13.5751, 23.5542, 40.7241, 70.5309, 122.164,
              211.609, 366.646)
grid_orientation = (3.53547, 2.75007, 3.53569, 2.74962, 3.92596, 2.35647,
                    0.786083, 0.00140674, 1.57147)

data_servers = ('http://hesperia.gsfc.nasa.gov/hessidata/',
                'http://hessi.ssl.berkeley.edu/hessidata/',
                'http://soleil.i4ds.ch/hessidata/')

lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                 'navy', 'orange', 'green')


def get_obssumm_dbase_file(time_range):
    """
    Download the RHESSI observing summary database file. This file lists the
    name of observing summary files for specific time ranges.

    Parameters
    ----------
    time_range : `str`, `sunpy.time.TimeRange`
        A `~sunpy.time.TimeRange` or `~sunpy.time.TimeRange` compatible string.

    Returns
    -------
    value : `tuple`
        Return a `tuple` (filename, headers) where filename is the local file
        name under which the object can be found, and headers is
        whatever the info() method of the object returned by urlopen.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))   # doctest: +SKIP

    References
    ----------
    | http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing Summary Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """

    #    http://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt

    _time_range = TimeRange(time_range)
    data_location = 'dbase/'

    url_root = data_servers[0] + data_location
    url = url_root + _time_range.start.strftime("hsi_obssumm_filedb_%Y%m.txt")

    f = urllib.request.urlretrieve(url)

    return f


def parse_obssumm_dbase_file(filename):
    """
    Parse the RHESSI observing summary database file. This file lists the
    name of observing summary files for specific time ranges along with other
    info

    Parameters
    ----------
    filename : `str`
        The filename of the obssumm dbase file.

    Returns
    -------
    out : `dict`
        Return a `dict` containing the parsed data in the dbase file.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> f = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))   # doctest: +SKIP
    >>> rhessi.parse_obssumm_dbase_file(f[0])   # doctest: +SKIP

    References
    ----------
    | http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing Summary Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    with open(filename, "rt") as fd:
        reader = csv.reader(fd, delimiter=' ', skipinitialspace=True)
        headerline = next(reader)
        headerline = next(reader)
        headerline = next(reader)
        headerline = next(reader)

        obssumm_filename = []
        orbit_start = []
        orbit_end = []
        start_time = []
        end_time = []
        status_flag = []
        number_of_packets = []

        for row in reader:
            obssumm_filename.append(row[0])
            orbit_start.append(int(row[1]))
            orbit_end.append(int(row[2]))
            start_time.append(datetime.strptime(row[3], '%d-%b-%y'))
            end_time.append(datetime.strptime(row[5], '%d-%b-%y'))
            status_flag.append(int(row[7]))
            number_of_packets.append(int(row[8]))

        return {
            headerline[0].lower(): obssumm_filename,
            headerline[1].lower(): orbit_start,
            headerline[2].lower(): orbit_end,
            headerline[3].lower(): start_time,
            headerline[4].lower(): end_time,
            headerline[5].lower(): status_flag,
            headerline[6].lower(): number_of_packets
        }


def get_obssum_filename(time_range):
    """
    Download the RHESSI observing summary data from one of the RHESSI
    servers, parses it, and returns the name of the obssumm file relevant for
    the time range

    Parameters
    ----------
    time_range : str, TimeRange
        A TimeRange or time range compatible string

    Returns
    -------
    out : string
        Returns the filename of the observation summary file

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_filename(('2011/04/04', '2011/04/05'))   # doctest: +SKIP

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    # need to download and inspect the dbase file to determine the filename
    # for the observing summary data
    f = get_obssumm_dbase_file(time_range)
    data_location = 'metadata/catalog/'

    result = parse_obssumm_dbase_file(f[0])
    _time_range = TimeRange(time_range)

    index_number = _time_range.start.day - 1

    return data_servers[0] + data_location + result.get('filename')[index_number] + 's'


def get_obssumm_file(time_range):
    """
    Download the RHESSI observing summary data from one of the RHESSI
    servers.

    Parameters
    ----------
    time_range : `str`, `sunpy.time.TimeRange`
        A TimeRange or time range compatible string

    Returns
    -------
    out : tuple
        Return a tuple (filename, headers) where filename is the local file
        name under which the object can be found, and headers is
        whatever the info() method of the object returned by urlopen.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))   # doctest: +SKIP

    .. note::
        This API is currently limited to providing data from whole days only.

    """

    time_range = TimeRange(time_range)
    data_location = 'metadata/catalog/'

    # TODO need to check which is the closest servers
    url_root = data_servers[0] + data_location

    url = url_root + get_obssum_filename(time_range)

    print('Downloading file: ' + url)
    f = urllib.request.urlretrieve(url)

    return f


def parse_obssumm_file(filename):
    """
    Parse a RHESSI observation summary file.

    Parameters
    ----------
    filename : str
        The filename of a RHESSI fits file.

    Returns
    -------
    out : `dict`
        Returns a dictionary.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> f = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))   # doctest: +SKIP
    >>> data = rhessi.parse_obssumm_file(f[0])   # doctest: +SKIP

    """

    afits = fits.open(filename)
    header = afits[0].header

    reference_time_ut = parse_time(afits[5].data.field('UT_REF')[0])
    time_interval_sec = afits[5].data.field('TIME_INTV')[0]
    # label_unit = fits[5].data.field('DIM1_UNIT')[0]
    # labels = fits[5].data.field('DIM1_IDS')
    labels = ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV',
              '50 - 100 keV', '100 - 300 keV', '300 - 800 keV', '800 - 7000 keV',
              '7000 - 20000 keV']

    # the data stored in the fits file are "compressed" countrates stored as one byte
    compressed_countrate = np.array(afits[6].data.field('countrate'))

    countrate = uncompress_countrate(compressed_countrate)
    dim = np.array(countrate[:,0]).size

    time_array = [reference_time_ut + timedelta(0,time_interval_sec * a) for a in np.arange(dim)]

    #TODO generate the labels for the dict automatically from labels
    data = {'time': time_array, 'data': countrate, 'labels': labels}

    return header, data

def uncompress_countrate(compressed_countrate):
    """Convert the compressed count rate inside of observing summary file from
    a compressed byte to a true count rate

    Parameters
    ----------
    compressed_countrate : byte array
        A compressed count rate returned from an observing summary file.

    References
    ----------
    Hsi_obs_summ_decompress.pro `<http://hesperia.gsfc.nasa.gov/ssw/hessi/idl/qlook_archive/hsi_obs_summ_decompress.pro>`_
    """
    ll = np.arange(0, 16, 1)
    lkup = np.zeros(256, dtype='int')
    sum = 0
    for i in range(0, 16):
        lkup[16 * i:16 * (i + 1)] = ll * 2 ** i + sum
        if i < 15:
            sum = lkup[16 * (i + 1) - 1] + 2 ** i
    return lkup[compressed_countrate]

def hsi_linecolors():
    """Define discrete colors to use for RHESSI plots

    Parameters
    ----------
    None

    Returns
    -------
    tuple : matplotliblib color list

    References
    ----------
    hsi_linecolors.pro `<http://hesperia.gsfc.nasa.gov/ssw/hessi/idl/gen/hsi_linecolors.pro`_
    """
    return ('black', 'magenta', 'lime', 'cyan', 'y', 'red', 'blue', 'orange', 'olive')

def _backproject(calibrated_event_list, detector=8, pixel_size=(1., 1.),
                 image_dim=(64, 64)):
    """
    Given a stacked calibrated event list fits file create a back
    projection image for an individual detectors. This function is used by
    backprojection.

    Parameters
    ----------
    calibrated_event_list : string
        filename of a RHESSI calibrated event list
    detector : int
        the detector number
    pixel_size : 2-tuple
        the size of the pixels in arcseconds. Default is (1,1).
    image_dim : 2-tuple
        the size of the output image in number of pixels

    Returns
    -------
    out : ndarray
        Return a backprojection image.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi

    """
    afits = fits.open(calibrated_event_list)

    #info_parameters = fits[2]
    #detector_efficiency = info_parameters.data.field('cbe_det_eff$$REL')

    afits = fits.open(calibrated_event_list)

    fits_detector_index = detector + 2
    detector_index = detector - 1
    grid_angle = np.pi/2. - grid_orientation[detector_index]
    harm_ang_pitch = grid_pitch[detector_index]/1

    phase_map_center = afits[fits_detector_index].data.field('phase_map_ctr')
    this_roll_angle = afits[fits_detector_index].data.field('roll_angle')
    modamp = afits[fits_detector_index].data.field('modamp')
    grid_transmission = afits[fits_detector_index].data.field('gridtran')
    count = afits[fits_detector_index].data.field('count')

    tempa = (np.arange(image_dim[0]*image_dim[1]) %  image_dim[0]) - (image_dim[0]-1)/2.
    tempb = tempa.reshape(image_dim[0],image_dim[1]).transpose().reshape(image_dim[0]*image_dim[1])

    pixel = np.array(list(zip(tempa,tempb)))*pixel_size[0]
    phase_pixel = (2*np.pi/harm_ang_pitch)* ( np.outer(pixel[:,0], np.cos(this_roll_angle - grid_angle)) -
                                              np.outer(pixel[:,1], np.sin(this_roll_angle - grid_angle))) + phase_map_center
    phase_modulation = np.cos(phase_pixel)
    gridmod = modamp * grid_transmission
    probability_of_transmission = gridmod*phase_modulation + grid_transmission
    bproj_image = np.inner(probability_of_transmission, count).reshape(image_dim)

    return bproj_image

def backprojection(calibrated_event_list, pixel_size=(1., 1.) * u.arcsec,
                   image_dim=(64, 64) * u.pix):
    """
    Given a stacked calibrated event list fits file create a back
    projection image.

    .. warning:: The image is not in the right orientation!

    Parameters
    ----------
    calibrated_event_list : string
        filename of a RHESSI calibrated event list
    pixel_size : `~astropy.units.Quantity` instance
        the size of the pixels in arcseconds. Default is (1,1).
    image_dim : `~astropy.units.Quantity` instance
        the size of the output image in number of pixels

    Returns
    -------
    out : RHESSImap
        Return a backprojection map.

    Examples
    --------
    >>> import sunpy.data
    >>> import sunpy.data.sample
    >>> import sunpy.instr.rhessi as rhessi
    >>> sunpy.data.download_sample_data(overwrite=False)   # doctest: +SKIP
    >>> map = rhessi.backprojection(sunpy.data.sample.RHESSI_EVENT_LIST)   # doctest: +SKIP
    >>> map.peek()   # doctest: +SKIP

    """
    if not isinstance(pixel_size, u.Quantity):
        raise ValueError("Must be astropy Quantity in arcseconds")
    try:
        pixel_size = pixel_size.to(u.arcsec)
    except:
        raise ValueError("'{0}' is not a valid pixel_size unit".format(pixel_size.unit))
    if not (isinstance(image_dim, u.Quantity) and image_dim.unit == 'pix'):
        raise ValueError("Must be astropy Quantity in pixels")

    try:
        import sunpy.data.sample
    except ImportError:
        import sunpy.data
        sunpy.data.download_sample()
    # This may need to be moved up to data from sample
    calibrated_event_list = sunpy.data.sample.RHESSI_EVENT_LIST

    afits = fits.open(calibrated_event_list)
    info_parameters = afits[2]
    xyoffset = info_parameters.data.field('USED_XYOFFSET')[0]
    time_range = TimeRange(info_parameters.data.field('ABSOLUTE_TIME_RANGE')[0])

    image = np.zeros(image_dim.value)

    # find out what detectors were used
    det_index_mask = afits[1].data.field('det_index_mask')[0]
    detector_list = (np.arange(9)+1) * np.array(det_index_mask)
    for detector in detector_list:
        if detector > 0:
            image = image + _backproject(calibrated_event_list, detector=detector, pixel_size=pixel_size.value
										 , image_dim=image_dim.value)

    dict_header = {
        "DATE-OBS": time_range.center().strftime("%Y-%m-%d %H:%M:%S"),
        "CDELT1": pixel_size[0],
        "NAXIS1": image_dim[0],
        "CRVAL1": xyoffset[0],
        "CRPIX1": image_dim[0].value/2 + 0.5,
        "CUNIT1": "arcsec",
        "CTYPE1": "HPLN-TAN",
        "CDELT2": pixel_size[1],
        "NAXIS2": image_dim[1],
        "CRVAL2": xyoffset[1],
        "CRPIX2": image_dim[0].value/2 + 0.5,
        "CUNIT2": "arcsec",
        "CTYPE2": "HPLT-TAN",
        "HGLT_OBS": 0,
        "HGLN_OBS": 0,
        "RSUN_OBS": solar_semidiameter_angular_size(time_range.center()).value,
        "RSUN_REF": sunpy.sun.constants.radius.value,
        "DSUN_OBS": sunearth_distance(time_range.center()) * sunpy.sun.constants.au.value
    }

    header = sunpy.map.MapMeta(dict_header)
    result_map = sunpy.map.Map(image, header)

    return result_map
