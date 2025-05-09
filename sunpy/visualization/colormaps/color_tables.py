"""
This module provides dictionaries for generating
`~matplotlib.colors.LinearSegmentedColormap`, and a dictionary of these
dictionaries.
"""
import pathlib

import matplotlib.colors as colors
import numpy as np

import astropy.units as u

__all__ = [
    'aia_color_table', 'sswidl_lasco_color_table', 'eit_color_table',
    'sxt_color_table', 'xrt_color_table', 'trace_color_table',
    'sot_color_table', 'hmi_mag_color_table', 'suvi_color_table',
    'rhessi_color_table', 'std_gamma_2', 'euvi_color_table', 'solohri_lya1216_color_table',
     'suit_color_table', 'punch_color_table',
]


CMAP_DATA_DIR = pathlib.Path(__file__).parent.absolute() / 'data'


def create_cdict(r, g, b):
    """
    Create the color tuples in the correct format.
    """
    i = np.linspace(0, 1, r.size)
    cdict = {name: list(zip(i, el / 255.0, el / 255.0))
             for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
    return cdict


def _cmap_from_rgb(r, g, b, name):
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(name, cdict)


def cmap_from_rgb_file(name, fname):
    """
    Create a colormap from a RGB .csv file.

    The .csv file must have 3  equal-length columns of integer data, with values
    between 0 and 255, which are the red, green, and blue values for the colormap.

    Parameters
    ----------
    name : str
        Name of the colormap.
    fname : str
        Filename of data file. Relative to the sunpy colormap data directory.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    data = np.loadtxt(CMAP_DATA_DIR / fname, delimiter=',')
    if data.shape[1] != 3:
        raise RuntimeError(f'RGB data files must have 3 columns (got {data.shape[1]})')
    return _cmap_from_rgb(data[:, 0], data[:, 1], data[:, 2], name)


def get_idl3():
    # The following values describe color table 3 for IDL (Red Temperature)
    return np.loadtxt(CMAP_DATA_DIR / 'idl_3.csv', delimiter=',')


def solohri_lya1216_color_table():
    solohri_lya1216 = get_idl3()
    solohri_lya1216[:, 2] = solohri_lya1216[:, 0] * np.linspace(0, 1, 256)
    return _cmap_from_rgb(*solohri_lya1216.T, 'SolO EUI HRI Lyman Alpha')


def create_aia_wave_dict():
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]

    c0 = np.arange(256, dtype='f')
    c1 = (np.sqrt(c0) * np.sqrt(255.0)).astype('f')
    c2 = (np.arange(256)**2 / 255.0).astype('f')
    c3 = ((c1 + c2 / 2.0) * 255.0 / (c1.max() + c2.max() / 2.0)).astype('f')

    aia_wave_dict = {
        1600*u.angstrom: (c3, c3, c2),
        1700*u.angstrom: (c1, c0, c0),
        4500*u.angstrom: (c0, c0, b0 / 2.0),
        94*u.angstrom: (c2, c3, c0),
        131*u.angstrom: (g0, r0, r0),
        171*u.angstrom: (r0, c0, b0),
        193*u.angstrom: (c1, c0, c2),
        211*u.angstrom: (c1, c0, c3),
        304*u.angstrom: (r0, g0, b0),
        335*u.angstrom: (c2, c0, c1)
    }
    return aia_wave_dict


@u.quantity_input
def aia_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SDO AIA images.

    Based on aia_lct.pro part of SDO/AIA on SSWIDL written by Karel
    Schrijver (2010/04/12).

    Parameters
    ----------
    wavelength : `~astropy.units.quantity.Quantity`
        Wavelength for the desired AIA color table.
    """
    aia_wave_dict = create_aia_wave_dict()
    try:
        r, g, b = aia_wave_dict[wavelength]
    except KeyError:
        raise ValueError("Invalid AIA wavelength. Valid values are "
                         "1600,1700,4500,94,131,171,193,211,304,335.")

    return _cmap_from_rgb(r, g, b, f'SDO AIA {str(wavelength):s}')


@u.quantity_input
def eit_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SOHO EIT images.
    """
    # SOHO EIT Color tables
    # EIT 171 IDL Name EIT Dark Bot Blue
    # EIT 195 IDL Name EIT Dark Bot Green
    # EIT 284 IDL Name EIT Dark Bot Yellow
    # EIT 304 IDL Name EIT Dark Bot Red
    try:
        color = {171*u.angstrom: 'dark_blue', 195*u.angstrom: 'dark_green',
                 284*u.angstrom: 'yellow', 304*u.angstrom: 'dark_red',
                 }[wavelength]
    except KeyError:
        raise ValueError("Invalid EIT wavelength. Valid values are "
                         "171, 195, 284, 304.")

    return cmap_from_rgb_file(f'SOHO EIT {str(wavelength):s}', f'eit_{color}.csv')


def sswidl_lasco_color_table(number):
    """
    Returns one of the SSWIDL-defined color tables for SOHO LASCO images.

    This function is included to allow users to access the SSWIDL-
    defined LASCO color tables provided by SunPy. It is recommended to
    use the function 'lasco_color_table' to obtain color tables for use
    with LASCO data and Helioviewer JP2 images.
    """
    try:
        return cmap_from_rgb_file(f'SOHO LASCO C{number}', f'lasco_c{number}.csv')
    except OSError:
        raise ValueError("Invalid LASCO number. Valid values are 2, 3.")


# Translated from the JP2Gen IDL SXT code lct_yla_gold.pro. Might be better
# to explicitly copy the numbers from the IDL calculation. This is a little
# more compact.
sxt_gold_r = np.concatenate((np.linspace(0, 255, num=185,
                                         endpoint=False), 255 * np.ones(71)))
sxt_gold_g = 255 * (np.arange(256)**1.25) / (255.0**1.25)
sxt_gold_b = np.concatenate((np.zeros(185), 255.0 * np.arange(71) / 71.0))

grayscale = np.arange(256)
grayscale.setflags(write=False)


def sxt_color_table(sxt_filter):
    """
    Returns one of the fundamental color tables for Yokhoh SXT images.
    """
    try:
        r, g, b = {
            'al': (sxt_gold_r, sxt_gold_g, sxt_gold_b),
            'wh': (grayscale, grayscale, grayscale)
        }[sxt_filter]
    except KeyError:
        raise ValueError("Invalid SXT filter type number. Valid values are "
                         "'al', 'wh'.")
    return _cmap_from_rgb(r, g, b, f'Yohkoh SXT {sxt_filter.title():s}')


def xrt_color_table():
    """
    Returns the color table used for all Hinode XRT images.
    """
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]
    return _cmap_from_rgb(r0, g0, b0, 'Hinode XRT')


def cor_color_table(number):
    """
    Returns one of the fundamental color tables for STEREO coronagraph images.
    """
    # STEREO COR Color tables
    if number not in [1, 2]:
        raise ValueError("Invalid COR number. Valid values are " "1, 2.")

    return cmap_from_rgb_file(f'STEREO COR{number}', f'stereo_cor{number}.csv')


def trace_color_table(measurement):
    """
    Returns one of the standard color tables for TRACE JP2 files.
    """
    if measurement == 'WL':
        return cmap_from_rgb_file(f'TRACE {measurement}', 'grayscale.csv')

    try:
        return cmap_from_rgb_file(f'TRACE {measurement}', f'trace_{measurement}.csv')
    except OSError:
        raise ValueError(
            "Invalid TRACE filter waveband passed. Valid values are "
            "171, 195, 284, 1216, 1550, 1600, 1700, WL")


def sot_color_table(measurement):
    """
    Returns one of the standard color tables for SOT files (following osdc
    convention).

    The relations between observation and color have been defined in
    hinode.py
    """
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]
    try:
        r, g, b = {
            'intensity': (r0, g0, b0),
        }[measurement]
    except KeyError:
        raise ValueError(
            "Invalid (or not supported) SOT type. Valid values are: "
            "intensity")

    return _cmap_from_rgb(r, g, b, f'Hinode SOT {measurement:s}')


def iris_sji_color_table(measurement, aialike=False):
    """
    Return the standard color table for IRIS SJI files.
    """
    # base vectors for IRIS SJI color tables
    c0 = np.arange(0, 256)
    c1 = (np.sqrt(c0) * np.sqrt(255)).astype(np.uint8)
    c2 = (c0**2 / 255.).astype(np.uint8)
    c3 = ((c1 + c2 / 2.) * 255. / (np.max(c1) + np.max(c2) / 2.)).astype(
        np.uint8)
    c4 = np.zeros(256).astype(np.uint8)
    c4[50:256] = (1 / 165. * np.arange(0, 206)**2).astype(np.uint8)
    c5 = ((1 + c1 + c3.astype(np.uint)) / 2.).astype(np.uint8)

    rr = np.ones(256, dtype=np.uint8) * 255
    rr[0:176] = np.arange(0, 176) / 175. * 255.
    gg = np.zeros(256, dtype=np.uint8)
    gg[100:256] = np.arange(0, 156) / 155. * 255.
    bb = np.zeros(256, dtype=np.uint8)
    bb[150:256] = np.arange(0, 106) / 105. * 255.
    agg = np.zeros(256, dtype=np.uint8)
    agg[120:256] = np.arange(0, 136) / 135. * 255.
    abb = np.zeros(256, dtype=np.uint8)
    abb[190:256] = np.arange(0, 66) / 65. * 255.

    if aialike:
        color_table = {
            '1330': (c1, c0, c2),
            '1400': (rr, agg, abb),
            '2796': (rr, c0, abb),
            '2832': (c3, c3, c2),
        }
    else:
        color_table = {
            '1330': (rr, gg, bb),
            '1400': (c5, c2, c4),
            '2796': (c1, c3, c2),
            '2832': (c0, c0, c2),
        }

    color_table.update({
        '1600': (c1, c0, c0),
        '5000': (c1, c1, c0),
        'FUV': (rr, gg, bb),
        'NUV': (c1, c3, c2),
        'SJI_NUV': (c0, c0, c0)
    })

    try:
        r, g, b = color_table[measurement]
    except KeyError:
        raise ValueError("Invalid IRIS SJI waveband. Valid values are \n" +
                         str(list(color_table.keys())))

    return _cmap_from_rgb(r, g, b, f'IRIS SJI {measurement:s}')


def hmi_mag_color_table():
    """
    Returns an alternate HMI Magnetogram color table; from Stanford
    University/JSOC.

    This is used by default for the `~sunpy.map.sources.HMISynopticMap`.
    But by default, for `~sunpy.map.sources.HMIMap` a grayscale colormap is used.

    Examples
    --------

    .. plot::
        :include-source:
        :context: close-figs

        import matplotlib.pyplot as plt
        import astropy.units as u

        import sunpy.map
        from sunpy.data.sample import HMI_LOS_IMAGE

        smap = sunpy.map.Map(HMI_LOS_IMAGE)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection=smap)

        smap.plot(axes=ax, cmap="hmimag", norm=None, vmin=-1500.0, vmax=1500.0)

        plt.show()

    References
    ----------
    * `Stanford Colortable (pdf) <http://jsoc.stanford.edu/data/hmi/HMI_M.ColorTable.pdf>`__
    """
    return cmap_from_rgb_file('SDO HMI magnetogram', 'hmi_mag.csv')


def stereo_hi_color_table(camera):
    if camera not in [1, 2]:
        raise ValueError("Valid HI cameras are 1 and 2")
    return cmap_from_rgb_file(f'STEREO HI{camera}', f'hi{camera}.csv')


@u.quantity_input
def suvi_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SUVI images.
    SUVI uses AIA color tables.
    """
    aia_wave_dict = create_aia_wave_dict()
    try:
        if wavelength == 195*u.angstrom:
            r, g, b = aia_wave_dict[193*u.angstrom]
        elif wavelength == 284*u.angstrom:
            r, g, b = aia_wave_dict[335*u.angstrom]
        else:
            r, g, b = aia_wave_dict[wavelength]
    except KeyError:
        raise ValueError(
            "Invalid SUVI wavelength. Valid values are "
            "94, 131, 171, 195, 284, 304."
        )
    return _cmap_from_rgb(r, g, b, f'GOES-R SUVI {str(wavelength):s}')


def rhessi_color_table():
    return cmap_from_rgb_file("rhessi", "rhessi.csv")


def std_gamma_2():
    return cmap_from_rgb_file("std_gamma_2", "std_gamma_2.csv")


def euvi_color_table(wavelength: u.angstrom):
    try:
        return cmap_from_rgb_file(f'EUVI {str(wavelength)}', f'euvi_{int(wavelength.value)}.csv')
    except OSError:
        raise ValueError(
            "Invalid EUVI wavelength. Valid values are "
            "171, 195, 284, 304."
        )


def suit_color_table(band):
    try:
        return cmap_from_rgb_file(f'suit_{band.lower()}', f'suit_{band.lower()}.csv')
    except OSError:
        raise ValueError(
                "Invalid Band. Valid Values are: "
                "NB01, NB02, NB03, NB04, NB05, NB06, NB07, NB08, BB01, BB02, BB03"
                )


def punch_color_table():
    return cmap_from_rgb_file("PUNCH", "punch.csv")
