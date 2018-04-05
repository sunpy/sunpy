# -*- coding: utf-8 -*-
"""
    Provides programs to process and analyze RHESSI data.

    .. warning:: This module is in development.

"""
from __future__ import absolute_import, print_function

import csv
import posixpath
import re
import socket
import warnings
from datetime import datetime, timedelta

import numpy as np
from dateutil.relativedelta import relativedelta

from astropy import units as u

from sunpy.time import TimeRange, parse_time
from sunpy.sun.sun import solar_semidiameter_angular_size
from sunpy.coordinates import get_sunearth_distance
import sunpy.map
import sunpy.io

from sunpy.extern.six.moves import urllib
from sunpy.extern.six.moves.urllib.request import urlopen, urlretrieve
from sunpy.extern.six.moves.urllib.error import URLError


__all__ = ['get_obssumm_dbase_file', 'parse_obssumm_dbase_file',
           'get_obssum_filename', 'get_obssumm_file', 'parse_obssumm_file',
           'backprojection']

# Measured fixed grid parameters
grid_pitch = (4.52467, 7.85160, 13.5751, 23.5542, 40.7241, 70.5309, 122.164,
              211.609, 366.646)
grid_orientation = (3.53547, 2.75007, 3.53569, 2.74962, 3.92596, 2.35647,
                    0.786083, 0.00140674, 1.57147)

data_servers = ('https://hesperia.gsfc.nasa.gov/hessidata/',
                'http://hessi.ssl.berkeley.edu/hessidata/',
                'http://soleil.i4ds.ch/hessidata/')

lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                 'navy', 'orange', 'green')


def get_base_url():
    """
    Find the first mirror which is online
    """
    for server in data_servers:
        try:
            urlopen(server, timeout=1)
            return server
        except (URLError, socket.timeout):
            pass

    raise IOError('Unable to find an online HESSI server from {0}'.format(data_servers))


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
    >>> fname, hdrs = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))   # doctest: +REMOTE_DATA

    References
    ----------
    | https://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing%20Summary%20Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    _time_range = TimeRange(time_range)

    if _time_range.start < parse_time("2002/02/01"):
        raise ValueError("RHESSI summary files are not available for before 2002-02-01")

    _check_one_day(_time_range)

    url = posixpath.join(get_base_url(), 'dbase',
                         _time_range.start.strftime("hsi_obssumm_filedb_%Y%m.txt"))

    return urlretrieve(url)


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
    >>> fname, _ = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))   # doctest: +REMOTE_DATA
    >>> file_names = rhessi.parse_obssumm_dbase_file(fname)   # doctest: +REMOTE_DATA
    >>> file_names['filename'][::5]   # doctest: +REMOTE_DATA
    ['hsi_obssumm_20110401_043.fit', 'hsi_obssumm_20110406_041.fit', 'hsi_obssumm_20110411_024.fit', 'hsi_obssumm_20110416_016.fit', 'hsi_obssumm_20110421_025.fit', 'hsi_obssumm_20110426_022.fit']

    References
    ----------
    | https://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing%20Summary%20Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    # An example dbase file can be found at:
    # https://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt

    with open(filename) as fd:
        reader = csv.reader(fd, delimiter=' ', skipinitialspace=True)
        _ = next(reader)  # skip 'HESSI Filedb File:' row
        _ = next(reader)  # skip 'Created: ...' row
        _ = next(reader)  # skip 'Number of Files: ...' row
        column_names = next(reader)  # ['Filename', 'Orb_st', 'Orb_end',...]

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
            start_time.append(datetime.strptime(row[3], '%d-%b-%y'))  # skip time
            end_time.append(datetime.strptime(row[5], '%d-%b-%y'))  # skip time
            status_flag.append(int(row[7]))
            number_of_packets.append(int(row[8]))

        return {
            column_names[0].lower(): obssumm_filename,
            column_names[1].lower(): orbit_start,
            column_names[2].lower(): orbit_end,
            column_names[3].lower(): start_time,
            column_names[4].lower(): end_time,
            column_names[5].lower(): status_flag,
            column_names[6].lower(): number_of_packets
        }


def get_obssum_filename(time_range):
    """
    Download the RHESSI observing summary data from one of the RHESSI
    servers, parses it, and returns the name of the obssumm files relevant for
    the time range.

    Parameters
    ----------
    time_range : str, TimeRange
        A TimeRange or time range compatible string

    Returns
    -------
    out : list
        Returns the filenames of the observation summary file

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssum_filename(('2011/04/04', '2011/04/05'))   # doctest: +REMOTE_DATA
    ['https://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20110404_042.fits']

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    time_range = TimeRange(time_range)

    delta = relativedelta(time_range.end, time_range.start)
    if delta.years > 0 or delta.months > 0:
        raise ValueError("Rhessi search results can not be found for a"
                         " time range crossing multiple months.")


    # need to download and inspect the dbase file to determine the filename
    # for the observing summary data

    dbase_file_name, _ = get_obssumm_dbase_file(time_range)
    dbase_dat = parse_obssumm_dbase_file(dbase_file_name)

    index_number_start = time_range.start.day - 1
    # If end is 0 set it to 1 so we always have at least one record.
    index_number_end = time_range.end.day - 1 or index_number_start + 1

    filenames = dbase_dat.get('filename')[index_number_start:index_number_end]
    return [posixpath.join(get_base_url(), 'metadata', 'catalog', filename + 's')
            for filename in filenames]


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
    >>> fname, hdrs = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))   # doctest: +REMOTE_DATA

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    _check_one_day(TimeRange(time_range))

    filenames = get_obssum_filename(time_range)

    # As we only support providing data from one whole day, only get the first file
    return urlretrieve(filenames[0])


def parse_obssumm_file(filename):
    """
    Parse a RHESSI observation summary file.
    Note: this is for the Lightcurve datatype only, the TimSeries uses the
    parse_obssumm_hdulist(hdulist) method to enable implicit source detection.

    Parameters
    ----------
    filename : str
        The filename of a RHESSI fits file.

    Returns
    -------
    value : `tuple`
        Return a `tuple` (fits_header, data). Where fits_header is of type
        `~astropy.io.fits.header.Header` and data of type `dict`

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> fname, _ = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))   # doctest: +REMOTE_DATA
    >>> data = rhessi.parse_obssumm_file(fname)   # doctest: +REMOTE_DATA

    """

    afits = sunpy.io.read_file(filename)
    fits_header = afits[0].header

    reference_time_ut = parse_time(afits[5].data.field('UT_REF')[0])
    time_interval_sec = afits[5].data.field('TIME_INTV')[0]

    # The data stored in the FITS file are "compressed" countrates stored as
    # one byte
    compressed_countrate = np.array(afits[6].data.field('countrate'))

    countrate = uncompress_countrate(compressed_countrate)
    dim = np.array(countrate[:, 0]).size

    time_array = [reference_time_ut + timedelta(0, time_interval_sec * a) for a in np.arange(dim)]

    labels = _build_energy_bands(label=afits[5].data.field('DIM1_UNIT')[0],
                                 bands=afits[5].data.field('DIM1_IDS')[0])

    return fits_header, dict(time=time_array, data=countrate, labels=labels)


def parse_obssumm_hdulist(hdulist):
    """
    Parse a RHESSI observation summary file.

    Parameters
    ----------
    hdulist : list
        The HDU list from the fits file.

    Returns
    -------
    out : `dict`
        Returns a dictionary.

    """
    header = hdulist[0].header

    reference_time_ut = parse_time(hdulist[5].data.field('UT_REF')[0])
    time_interval_sec = hdulist[5].data.field('TIME_INTV')[0]
    # label_unit = fits[5].data.field('DIM1_UNIT')[0]
    # labels = fits[5].data.field('DIM1_IDS')
    labels = ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV',
              '50 - 100 keV', '100 - 300 keV', '300 - 800 keV',
              '800 - 7000 keV', '7000 - 20000 keV']

    # The data stored in the fits file are "compressed" countrates stored as
    # one byte
    compressed_countrate = np.array(hdulist[6].data.field('countrate'))

    countrate = uncompress_countrate(compressed_countrate)
    dim = np.array(countrate[:, 0]).size

    time_array = [reference_time_ut + timedelta(0, time_interval_sec * a) for a in np.arange(dim)]

    #  TODO generate the labels for the dict automatically from labels
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
    Hsi_obs_summ_decompress.pro `<https://hesperia.gsfc.nasa.gov/ssw/hessi/idl/qlook_archive/hsi_obs_summ_decompress.pro>`_
    """

    # Ensure uncompressed counts are between 0 and 255
    if (compressed_countrate.min() < 0) or (compressed_countrate.max() > 255):
        raise ValueError(
            'Exepected uncompressed counts {} to in range 0-255'.format(compressed_countrate))

    # TODO Must be a better way than creating entire lookup table on each call
    ll = np.arange(0, 16, 1)
    lkup = np.zeros(256, dtype='int')
    _sum = 0
    for i in range(0, 16):
        lkup[16 * i:16 * (i + 1)] = ll * 2 ** i + _sum
        if i < 15:
            _sum = lkup[16 * (i + 1) - 1] + 2 ** i

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
    hsi_linecolors.pro `<https://hesperia.gsfc.nasa.gov/ssw/hessi/idl/gen/hsi_linecolors.pro>`_
    """
    return ('black', 'magenta', 'lime', 'cyan', 'y', 'red', 'blue', 'orange',
            'olive')


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
    # info_parameters = fits[2]
    # detector_efficiency = info_parameters.data.field('cbe_det_eff$$REL')

    afits = sunpy.io.read_file(calibrated_event_list)

    fits_detector_index = detector + 2
    detector_index = detector - 1
    grid_angle = np.pi/2. - grid_orientation[detector_index]
    harm_ang_pitch = grid_pitch[detector_index]/1

    phase_map_center = afits[fits_detector_index].data.field('phase_map_ctr')
    this_roll_angle = afits[fits_detector_index].data.field('roll_angle')
    modamp = afits[fits_detector_index].data.field('modamp')
    grid_transmission = afits[fits_detector_index].data.field('gridtran')
    count = afits[fits_detector_index].data.field('count')

    tempa = (np.arange(image_dim[0] * image_dim[1]) % image_dim[0]) - (image_dim[0]-1)/2.
    tempb = tempa.reshape(image_dim[0], image_dim[1]).transpose().reshape(image_dim[0]*image_dim[1])

    pixel = np.array(list(zip(tempa, tempb)))*pixel_size[0]
    phase_pixel = (2 * np.pi/harm_ang_pitch) *\
                  (np.outer(pixel[:, 0], np.cos(this_roll_angle - grid_angle)) -
                   np.outer(pixel[:, 1], np.sin(this_roll_angle - grid_angle))) + phase_map_center
    phase_modulation = np.cos(phase_pixel)
    gridmod = modamp * grid_transmission
    probability_of_transmission = gridmod * phase_modulation + grid_transmission
    bproj_image = np.inner(probability_of_transmission, count).reshape(image_dim)

    return bproj_image


@u.quantity_input(pixel_size=u.arcsec, image_dim=u.pix)
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
    This example is broken.
    >>> import sunpy.data
    >>> import sunpy.data.sample  # doctest: +SKIP
    >>> import sunpy.instr.rhessi as rhessi
    >>> map = rhessi.backprojection(sunpy.data.sample.RHESSI_IMAGE)  # doctest: +SKIP
    >>> map.peek()   # doctest: +SKIP
    """
    pixel_size = pixel_size.to(u.arcsec)
    image_dim = np.array(image_dim.to(u.pix).value, dtype=int)

    afits = sunpy.io.read_file(calibrated_event_list)
    info_parameters = afits[2]
    xyoffset = info_parameters.data.field('USED_XYOFFSET')[0]
    time_range = TimeRange(info_parameters.data.field('ABSOLUTE_TIME_RANGE')[0])

    image = np.zeros(image_dim)

    # find out what detectors were used
    det_index_mask = afits[1].data.field('det_index_mask')[0]
    detector_list = (np.arange(9)+1) * np.array(det_index_mask)
    for detector in detector_list:
        if detector > 0:
            image = image + _backproject(calibrated_event_list, detector=detector,
                                         pixel_size=pixel_size.value, image_dim=image_dim)

    dict_header = {
        "DATE-OBS": time_range.center.strftime("%Y-%m-%d %H:%M:%S"),
        "CDELT1": pixel_size[0],
        "NAXIS1": image_dim[0],
        "CRVAL1": xyoffset[0],
        "CRPIX1": image_dim[0]/2 + 0.5,
        "CUNIT1": "arcsec",
        "CTYPE1": "HPLN-TAN",
        "CDELT2": pixel_size[1],
        "NAXIS2": image_dim[1],
        "CRVAL2": xyoffset[1],
        "CRPIX2": image_dim[0]/2 + 0.5,
        "CUNIT2": "arcsec",
        "CTYPE2": "HPLT-TAN",
        "HGLT_OBS": 0,
        "HGLN_OBS": 0,
        "RSUN_OBS": solar_semidiameter_angular_size(time_range.center).value,
        "RSUN_REF": sunpy.sun.constants.radius.value,
        "DSUN_OBS": get_sunearth_distance(time_range.center).value * sunpy.sun.constants.au.value
    }

    result_map = sunpy.map.Map(image, dict_header)

    return result_map


def _build_energy_bands(label, bands):
    """
    Parameters
    ----------
    label: `str`
    bands: `list` of `str`
    Returns
    -------
    bands_with_units: `list` of `str`
        Each `str` item is an energy band and its unit
    Example
    -------
    >>> from sunpy.instr.rhessi import _build_energy_bands
    >>> _build_energy_bands('Energy bands (keV)', ['3 - 6', '6 - 12', '12 - 25'])
    ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV']
    """

    unit_pattern = re.compile(r'^.+\((?P<UNIT>\w+)\)$')

    matched = unit_pattern.match(label)

    if matched is None:
        raise ValueError("Unable to find energy unit in '{0}' "
                         "using REGEX '{1}'".format(label, unit_pattern.pattern))

    unit = matched.group('UNIT').strip()

    return ['{energy_band} {unit}'.format(energy_band=band, unit=unit) for band in bands]


def _check_one_day(time_range):
    """
    Currently only support TimeRanges of a maximum of one day.
    Issue a visible warning if `time_range` is greater than this
    Parameters
    ----------
    time_range : `sunpy.time.TimeRange`
    """
    if time_range.days > 1 * u.day:
        warnings.warn('Currently only support providing data from one whole day. Only data for {0} '
                      'will be returned'.format(time_range.start.strftime("%Y-%m-%d")), UserWarning,
                      stacklevel=2)
