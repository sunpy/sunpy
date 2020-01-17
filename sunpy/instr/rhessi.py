"""
This module provides processing routines programs to process and analyze RHESSI
data.
"""

import re
import csv

from dateutil.relativedelta import relativedelta

import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time, TimeDelta

import sunpy.io
from sunpy.coordinates import sun
from sunpy.time import TimeRange, parse_time

__all__ = [
    'parse_observing_summary_hdulist', 'backprojection', 'parse_observing_summary_dbase_file',
    'get_flare_list', 'read_flare_list_file', 'print_flare_list', 'convert_flag_dict'
]


# Measured fixed grid parameters
grid_pitch = (4.52467, 7.85160, 13.5751, 23.5542, 40.7241, 70.5309, 122.164,
              211.609, 366.646)
grid_orientation = (3.53547, 2.75007, 3.53569, 2.74962, 3.92596, 2.35647,
                    0.786083, 0.00140674, 1.57147)

lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                 'navy', 'orange', 'green')

KNOWN_FLARE_LIST_SOURCES = {
    "NASA": "https://hesperia.gsfc.nasa.gov/hessidata/dbase/",
    "Berkeley": "http://hessi.ssl.berkeley.edu/hessidata/dbase/",
    "i4DS": "http://soleil.i4ds.ch/hessidata/dbase/",
}

FLAGID2FLAG = {
    'SAA_AT_START': 'SS',
    'SAA_AT_END': 'SE',
    'SAA_DURING_FLARE': 'SD',
    'ECLIPSE_AT_START': 'ES',
    'ECLIPSE_AT_END': 'EE',
    'ECLIPSE_DURING_FLARE': 'ED',
    'FLARE_AT_SOF': 'FS',
    'FLARE_AT_EOF': 'FE',
    'NON_SOLAR': 'NS',
    'FAST_RATE_MODE': 'FR',
    'FRONT_DECIMATION': 'DF',
    'ATT_STATE_AT_PEAK': '',  # 0-3 -> changes flag an to An
    'DATA_GAP_AT_START': 'GS',
    'DATA_GAP_AT_END': 'GE',
    'DATA_GAP_DURING_FLARE': 'GD',
    'PARTICLE_EVENT': 'PE',
    'DATA_QUALITY': 'Qn',  # 0-11
    'POSITION_QUALITY': 'P1',
    'ATTEN_0': 'a0',
    'ATTEN_1': 'a1',
    'ATTEN_2': 'a2',
    'ATTEN_3': 'a3',
    'REAR_DECIMATION': 'DR',
    'MAGNETIC_REGION': 'MR',
    'IMAGE_STATUS': '',
    'SPECTRUM_STATUS': '',
    'SOLAR_UNCONFIRMED': 'PS',
    'SOLAR': ''
}


def parse_observing_summary_dbase_file(filename):
    """
    Parse the RHESSI observing summary database file.

    This file lists the name of observing summary files
    for specific time ranges along with other info.

    .. note::
        This API is currently limited to providing data from whole days only.

    Parameters
    ----------
    filename : `str`
        The filename of the obssumm dbase file.

    Returns
    -------
    `dict`
        Return a `dict` containing the parsed data in the dbase file.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.parse_observing_summary_dbase_file(fname)   # doctest: +SKIP

    References
    ----------
    https://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing%20Summary%20Data
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
            start_time.append(Time.strptime(row[3], '%d-%b-%y'))  # skip time
            end_time.append(Time.strptime(row[5], '%d-%b-%y'))  # skip time
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


def parse_observing_summary_hdulist(hdulist):
    """
    Parse a RHESSI observation summary file.

    Parameters
    ----------
    hdulist : `list`
        The HDU list from the fits file.

    Returns
    -------
    out : `dict`
        Returns a dictionary.
    """
    header = hdulist[0].header

    reference_time_ut = parse_time(hdulist[5].data.field('UT_REF')[0],
                                   format='utime')
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

    time_array = parse_time(reference_time_ut) + \
        TimeDelta(time_interval_sec * np.arange(dim) * u.second)

    #  TODO generate the labels for the dict automatically from labels
    data = {'time': time_array, 'data': countrate, 'labels': labels}

    return header, data


def uncompress_countrate(compressed_countrate):
    """
    Convert the compressed count rate inside of observing summary file from a
    compressed byte to a true count rate.

    Parameters
    ----------
    compressed_countrate : `byte` array
        A compressed count rate returned from an observing summary file.

    References
    ----------
    `Hsi_obs_summ_decompress.pro <https://hesperia.gsfc.nasa.gov/ssw/hessi/idl/qlook_archive/hsi_obs_summ_decompress.pro>`_
    """

    # Ensure uncompressed counts are between 0 and 255
    if (compressed_countrate.min() < 0) or (compressed_countrate.max() > 255):
        raise ValueError(
            f'Exepected uncompressed counts {compressed_countrate} to in range 0-255')

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
    """
    Define discrete colors to use for RHESSI plots.

    Returns
    -------
    `tuple` :
         A tuple of names of colours.

    References
    ----------
    `hsi_linecolors.pro <https://hesperia.gsfc.nasa.gov/ssw/hessi/idl/gen/hsi_linecolors.pro>`__
    """
    return ('black', 'magenta', 'lime', 'cyan', 'y', 'red', 'blue', 'orange',
            'olive')


def _backproject(calibrated_event_list, detector=8, pixel_size=(1., 1.),
                 image_dim=(64, 64)):
    """
    Given a stacked calibrated event list fits file create a back projection
    image for an individual detectors.

    Parameters
    ----------
    calibrated_event_list : `str`
        Filename of a RHESSI calibrated event list.
    detector : `int`, optional
        The detector number.
    pixel_size : `tuple`, optional
        A length 2 tuple with the size of the pixels in arcseconds.
        Defaults to  ``(1, 1)``.
    image_dim : `tuple`, optional
        A length 2 tuple with the size of the output image in number of pixels.
        Defaults to ``(64, 64)``.

    Returns
    -------
    `numpy.ndarray`
        A backprojection image.
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


@u.quantity_input
def backprojection(calibrated_event_list, pixel_size: u.arcsec = (1., 1.) * u.arcsec,
                   image_dim: u.pix = (64, 64) * u.pix):
    """
    Given a stacked calibrated event list fits file create a back projection
    image.

    .. warning::

        The image will not be in the right orientation.

    Parameters
    ----------
    calibrated_event_list : `str`
        Filename of a RHESSI calibrated event list.
    pixel_size : `tuple`, optional
        A length 2 tuple with the size of the pixels in arcsecond
        `~astropy.units.Quantity`. Defaults to  ``(1, 1) * u.arcsec``.
    image_dim : `tuple`, optional
        A length 2 tuple with the size of the output image in number of pixel
        `~astropy.units.Quantity` Defaults to ``(64, 64) * u.pix``.

    Returns
    -------
    `sunpy.map.sources.RHESSImap`
        A backprojection map.
    """
    # import sunpy.map in here so that net and timeseries don't end up importing map
    import sunpy.map

    pixel_size = pixel_size.to(u.arcsec)
    image_dim = np.array(image_dim.to(u.pix).value, dtype=int)

    afits = sunpy.io.read_file(calibrated_event_list)
    info_parameters = afits[2]
    xyoffset = info_parameters.data.field('USED_XYOFFSET')[0]
    time_range = TimeRange(info_parameters.data.field('ABSOLUTE_TIME_RANGE')[0], format='utime')

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
        "RSUN_OBS": sun.angular_radius(time_range.center).value,
        "RSUN_REF": sunpy.sun.constants.radius.value,
        "DSUN_OBS": sun.earth_distance(time_range.center).value * sunpy.sun.constants.au.value
    }

    result_map = sunpy.map.Map(image, dict_header)

    return result_map


def _build_energy_bands(label, bands):
    """
    Creates a list of strings with the correct formatting for axis labels.

    Parameters
    ----------
    label: `str`
        The ``label`` to use as a basis.
    bands: `list` of `str`
        The bands to append to the ``label``.

    Returns
    -------
    `list` of `str`
        Each string is an energy band and its unit.

    Example
    -------
    >>> from sunpy.instr.rhessi import _build_energy_bands
    >>> _build_energy_bands('Energy bands (keV)', ['3 - 6', '6 - 12', '12 - 25'])
    ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV']
    """

    unit_pattern = re.compile(r'^.+\((?P<UNIT>\w+)\)$')

    matched = unit_pattern.match(label)

    if matched is None:
        raise ValueError("Unable to find energy unit in '{}' "
                         "using REGEX '{}'".format(label, unit_pattern.pattern))

    unit = matched.group('UNIT').strip()

    return [f'{band} {unit}' for band in bands]


def get_flare_list(start: str,
                   end: str,
                   source: str='NASA',
                   file_format: str="hessi_flare_list_%Y%m.fits",
                   inc=relativedelta(months=+1)):
    """
    Read and combine RHESSI flare lists from .fits files as specified with further parameters
    Supported date formats are the same as in ``sunpy.time``.

    Parameters
    ----------
    start : `str`
        Start date of period within which flares should be loaded.
    end : `str`
        End date of period within which flares should be loaded.
    source : `str`, optional
        Source from where .fits files should be loaded. Can be an URL or a local folder.
        Known sources are "NASA", "Berkeley" and "i4DS".
        Defaults to ``"NASA"``.
    file_format : `str`, optional
        Specifies the naming convention of the files available in the source folder.
        The changing parts should be given as format string as in ``strftime``.
        Defaults to ``"hessi_flare_list_%Y%m.fits"``.
    inc : `timedelta`, `relativedelta`, optional
        Specifies by how much time the single files are separated.
        Defaults to ``relativedelta(months=+1)``.

    Returns
    -------
    ``pandas.DataFrame``
        out : ``pandas.DataFrame`` containing the flares within the given time constraints

    Examples
    --------
    >>> from sunpy.instr.rhessi import get_flare_list
    >>> fl = get_flare_list("2018-01-06 16:32:57", "2018-01-22 02:43:27")  # doctest: +REMOTE_DATA
    """

    start_dt = parse_time(start).to_datetime()
    end_dt = parse_time(end).to_datetime()
    format_str = file_format[file_format.index("%"):file_format.rindex("%") + 2]
    cur_format = start_dt.strftime(format_str)
    end_format = end_dt.strftime(format_str)

    if source in KNOWN_FLARE_LIST_SOURCES:
        source = KNOWN_FLARE_LIST_SOURCES[source]

    cur_dt = start_dt
    result = pd.DataFrame()
    while cur_format <= end_format:
        file = file_format.replace(format_str, cur_format)
        result = result.append(read_flare_list_file(source + file), ignore_index=True)
        cur_dt = cur_dt + inc
        cur_format = cur_dt.strftime(format_str)

    # filter results for more detailed time constraints (if applicable)
    if len(end) <= 12:  # formats that do specify a time are at least 14 chars
        end_dt += relativedelta(days=+1)  # add day if end date was specified without time
    else:
        end_dt += relativedelta(microsecond=+1)  # add 1ms so "smaller" operator works as intended

    result = result[result['END_TIME'] >= start_dt]
    result = result[result['START_TIME'] < end_dt]
    return result


def read_flare_list_file(file):
    """
    Read RHESSI flare list .fits file into ``pandas.DataFrame``. TIME values are parsed with
    format 'utime', which is the same as Unix timestamp but starts 9 years later.
    FLAGS are assigned their respective label (FLAG ID) and returned as `dict`.

    Parameters
    ----------
    file : `str`
        The URL or local filename of the hessi flare list.

    Returns
    -------
    ``pandas.DataFrame``
        out : ``pandas.DataFrame`` containing the parsed flares.

    Examples
    --------
    >>> from sunpy.instr.rhessi import read_flare_list_file
    >>> url = "https://hesperia.gsfc.nasa.gov/hessidata/dbase/hessi_flare_list_201802.fits"
    >>> fl = read_flare_list_file(url)  # doctest: +REMOTE_DATA

    References
    ----------
    https://hesperia.gsfc.nasa.gov/rhessi3/data-access/rhessi-data/flare-list/index.html
    """
    try:
        fits = sunpy.io.fits.read(file)
    except Exception:
        raise RuntimeError("couldn't load file " + file)

    results = []
    for row in fits[3].data:
        result_row = {}
        for k in fits[3].data.columns.names:
            if k.endswith('_TIME'):
                result_row[k] = parse_time(row[k], format="utime")
                result_row[k].format = "datetime"  # for human readable display inside the DF
            elif k == 'FLAGS':
                flags = {}
                for i, fid in zip(row[k], fits[2].data['FLAG_IDS'][0]):
                    flags[fid] = i
                result_row[k] = flags
            else:
                result_row[k] = row[k]
        results.append(result_row)

    return pd.DataFrame(results)


def print_flare_list(data_frame):
    """
    Convert and print flares similar to the available .txt flare lists.

    Parameters
    ----------
    data_frame : ``pandas.DataFrame``
        DataFrame containing the flares to print.
    """

    for idx, row in data_frame.iterrows():
        print(
            "{id:9} {st} {pt} {et} {dur:5} {p:6} {n:9} {e:>11} {x:5} {y:5} {r:6} {ar:4}  {f}"
            .format(
                id=row['ID_NUMBER'],
                st=row['START_TIME'].strftime('%e-%b-%Y %H:%M:%S'),
                pt=row['PEAK_TIME'].strftime('%H:%M:%S'),
                et=row['END_TIME'].strftime('%H:%M:%S'),
                dur=int(round((row['END_TIME'] - row['START_TIME']).to_value("sec"))),
                p=int(row['PEAK_COUNTRATE']),
                n=int(row['TOTAL_COUNTS']),
                e=str(int(row['ENERGY_HI'][0])) + "-" + str(int(row['ENERGY_HI'][1])),
                x=int(round(row['POSITION'][0])),
                y=int(round(row['POSITION'][1])),
                r=int(round((row['POSITION'][0] ** 2 + row['POSITION'][1] ** 2) ** 0.5)),
                ar=row['ACTIVE_REGION'],
                f=" ".join(convert_flag_dict(row['FLAGS'])),
            )
        )


def convert_flag_dict(flags_dict):
    """
    Convert flag dictionary into array of short abbreviations (as in available .txt lists).

    Parameters
    ----------
    flags_dict : `dict`
        Dictionary containing the flags where key is the flag id and value the
        original value from the byte mask.

    Returns
    -------
    `list`
        out : `list` of 2-character strings, each representing a present flag
        (as in available .txt lists).

    """
    flags = []
    for k in flags_dict:
        # note: as of now the Q-flag is missing if the flag DATA_QUALITY is 0
        if flags_dict[k] > 0 and FLAGID2FLAG[k] != '':
            if k.startswith('ATTEN_') and int(k[-1:]) == flags_dict['ATT_STATE_AT_PEAK']:
                flags.append(FLAGID2FLAG[k].upper())
            else:
                flags.append(FLAGID2FLAG[k].replace("n", str(flags_dict[k])))
    flags.sort(key=str.lower)
    return flags
