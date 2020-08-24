"""
This module provides dictionaries for generating
`~matplotlib.colors.LinearSegmentedColormaps`, and a dictionary of these
dictionaries.
"""
import pathlib

import matplotlib.colors as colors
import numpy as np

import astropy.units as u

__all__ = [
    'aia_color_table', 'sswidl_lasco_color_table', 'eit_color_table',
    'sxt_color_table', 'xrt_color_table', 'trace_color_table',
    'sot_color_table', 'hmi_mag_color_table', 'suvi_color_table'
]

cmap_data_dir = pathlib.Path(__file__).parent.absolute() / 'data'

# The following values describe color table 3 for IDL (Red Temperature)
idl_3 = np.loadtxt(cmap_data_dir / 'idl_3.csv', delimiter=',')
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


@u.quantity_input
def aia_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SDO AIA images.

    Based on aia_lct.pro part of SDO/AIA on SSWIDL written by Karel
    Schrijver (2010/04/12).

    Parameters
    ----------
    wavelength : `~astropy.units.quantity`
        Wavelength for the desired AIA color table.
    """
    try:
        r, g, b = aia_wave_dict[wavelength]
    except KeyError:
        raise ValueError("Invalid AIA wavelength. Valid values are "
                         "1600,1700,4500,94,131,171,193,211,304,335.")

    # Now create the color dictionary in the correct format
    cdict = create_cdict(r, g, b)

    return colors.LinearSegmentedColormap(
        'SDO AIA {:s}'.format(str(wavelength)), cdict)


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

    return cmap_from_rgb_file('SOHO EIT {:s}'.format(str(wavelength)), f'eit_{color}.csv')


lasco_c2_r = np.concatenate((np.array([
    0., 1., 2., 5., 8., 11., 14., 17., 20., 23., 26., 28., 31., 34., 37., 42.,
    44., 47., 50., 55., 57., 60., 65., 68., 70., 75., 78., 82., 85., 88., 92.,
    95., 99., 102., 107., 110., 114., 117., 121., 124., 128., 133., 136., 140.,
    143., 147., 152., 155., 159., 163., 166., 170., 175., 178., 182., 186.,
    189., 194., 198., 201., 205., 210., 214., 217., 221., 226., 230., 233.,
    237., 241., 246., 250., 253.
]), 255 * np.ones(183)))

lasco_c2_g = np.concatenate(
    (np.zeros(52).astype('int'),
     np.array([
         1., 5., 11., 17., 20., 26., 32., 35., 41., 47., 52., 56., 62., 68.,
         73., 77., 83., 88., 94., 100., 103., 109., 115., 120., 126., 130.,
         136., 141., 147., 153., 158., 164., 168., 173., 179., 185., 190.,
         196., 202., 207., 213., 217., 222., 228., 234., 239., 245., 251.
     ]), 255 * np.ones(156)))

lasco_c2_b = np.concatenate(
    (np.zeros(78).astype('int'),
     np.array([
         0., 7., 19., 31., 43., 54., 66., 74., 86., 98., 109., 121., 133.,
         145., 156., 168., 176., 188., 200., 211., 223., 235., 247.
     ]), 255 * np.ones(156)))

lasco_c3_r = np.concatenate(
    (np.zeros(77).astype('int'),
     np.array([
         5., 13., 25., 33., 45., 53., 65., 73., 85., 94., 106., 114., 126.,
         134., 146., 154., 166., 175., 187., 195., 207., 215., 227., 235., 247.
     ]), 255 * np.ones(154)))

lasco_c3_g = np.concatenate(
    (np.zeros(39).astype('int'),
     np.array([
         4., 7., 12., 15., 20., 23., 28., 31., 36., 39., 44., 47., 52., 55.,
         60., 63., 68., 71., 76., 79., 84., 87., 92., 95., 100., 103., 108.,
         111., 116., 119., 124., 127., 132., 135., 140., 143., 148., 151.,
         156., 159., 164., 167., 172., 175., 180., 183., 188., 191., 196.,
         199., 204., 207., 212., 215., 220., 223., 228., 231., 236., 239.,
         244., 247., 252., 255., 255.
     ]), 255 * np.ones(154)))

lasco_c3_b = np.concatenate((np.array([
    0., 4., 6., 10., 13., 17., 20., 24., 27., 31., 33., 37., 40., 44., 47.,
    51., 54., 58., 61., 65., 67., 71., 74., 78., 81., 85., 88., 92., 94., 99.,
    101., 105., 108., 112., 115., 119., 122., 126., 128., 132., 135., 139.,
    142., 146., 149., 153., 155., 160., 162., 166., 169., 173., 176., 180.,
    183., 187., 189., 193., 196., 200., 203., 207., 210., 214., 217., 221.,
    223., 227., 230., 234., 237., 241., 244., 248., 250.
]), 255 * np.ones(181)))


def sswidl_lasco_color_table(number):
    """
    Returns one of the SSWIDL-defined color tables for SOHO LASCO images.

    This function is included to allow users to access the SSWIDL-
    defined LASCO color tables provided by SunPy. It is recommended to
    use the function 'lasco_color_table' to obtain color tables for use
    with LASCO data and Helioviewer JP2 images.
    """
    # SOHO LASCO Color tables
    # LASCO C2 white light IDL Name
    # LASCO C3 white light IDL Name
    try:
        r, g, b = {
            2: (lasco_c2_r, lasco_c2_g, lasco_c2_b),
            3: (lasco_c3_r, lasco_c3_g, lasco_c3_b),
        }[number]
    except KeyError:
        raise ValueError("Invalid LASCO number. Valid values are 2, 3.")

    # Now create the color dictionary in the correct format
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(
        'SOHO LASCO C{:s}'.format(str(number)), cdict)


# Translated from the JP2Gen IDL SXT code lct_yla_gold.pro.  Might be better
# to explicitly copy the numbers from the IDL calculation.  This is a little
# more compact.
sxt_gold_r = np.concatenate((np.linspace(0, 255, num=185,
                                         endpoint=False), 255 * np.ones(71)))
sxt_gold_g = 255 * (np.arange(256)**1.25) / (255.0**1.25)
sxt_gold_b = np.concatenate((np.zeros(185), 255.0 * np.arange(71) / 71.0))

grayscale = np.arange(256)


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

    # Now create the color dictionary in the correct format
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(
        'Yohkoh SXT {:s}'.format(sxt_filter.title()), cdict)


def xrt_color_table():
    """
    Returns the color table used for all Hinode XRT images.
    """
    # Now create the color dictionary in the correct format
    cdict = create_cdict(r0, g0, b0)
    return colors.LinearSegmentedColormap('Hinode XRT', cdict)


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
            "Invalid TRACE filter waveband passed.  Valid values are "
            "171, 195, 284, 1216, 1550, 1600, 1700, WL")


def sot_color_table(measurement):
    """
    Returns one of the standard color tables for SOT files (following osdc
    convention).

    The relations between observation and color have been defined in
    hinode.py
    """
    try:
        r, g, b = {
            'intensity': (r0, g0, b0),
        }[measurement]
    except KeyError:
        raise ValueError(
            "Invalid (or not supported) SOT type. Valid values are: "
            "intensity")

    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(f'Hinode SOT {measurement:s}', cdict)


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
        raise ValueError("Invalid IRIS SJI waveband.  Valid values are \n" +
                         str(list(color_table.keys())))

    # Now create the color dictionary in the correct format
    cdict = create_cdict(r, g, b)
    # Return the color table
    return colors.LinearSegmentedColormap(f'IRIS SJI {measurement:s}',
                                          cdict)


def hmi_mag_color_table():
    """
    Returns an alternate HMI Magnetogram color table; from Stanford
    University/JSOC.

    Examples
    --------
    >>> # Example usage for NRT data:
    >>> import sunpy.map
    >>> import matplotlib.pyplot as plt
    >>> hmi = sunpy.map.Map('fblos.fits')  # doctest: +SKIP
    >>> hmi.plot_settings['cmap'] = plt.get_cmap('hmimag')  # doctest: +SKIP
    >>> hmi.peek(vmin=-1500.0, vmax=1500.0)  # doctest: +SKIP

    >>> # OR (for a basic plot with pixel values on the axes)
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import sunpy.map
    >>> hmi = sunpy.map.Map('fblos.fits')  # doctest: +SKIP
    >>> plt.imshow(np.clip(hmi.data, -1500.0, 1500.0), cmap=plt.get_cmap('hmimag'), origin='lower')  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

    >>> # Example usage for science (Level 1.0) data:
    >>> import numpy as np
    >>> import sunpy.map
    >>> hmi = sunpy.map.Map('hmi.m_45s.2014.05.11_12_00_45_TAI.magnetogram.fits')  # doctest: +SKIP
    >>> hmir = hmi.rotate()  # doctest: +SKIP
    >>> hmir.plot_settings['cmap'] = plt.get_cmap('hmimag')  # doctest: +SKIP
    >>> hmir.peek(vmin=-1500.0, vmax=1500.0)  # doctest: +SKIP

    References
    ----------
    * `Stanford Colortable (pdf) <http://jsoc.stanford.edu/data/hmi/HMI_M.ColorTable.pdf>`_
    """
    return cmap_from_rgb_file('SDO HMI magnetogram', 'hmi_mag.csv')


def stereo_hi_color_table(camera):
    if camera in [1, 2]:
        return cmap_from_rgb_file(f'STEREO HI{camera}', f'hi{camera}.csv')
    else:
        raise ValueError("Valid HI cameras are 1 and 2")


def cmap_from_rgb_file(name, fname):
    """
    Create a colormap from a RGB .csv file.

    name : str
        Name of the colormap.
    fname : str
        Filename of data file. Relative to the sunpy colormap data directory.

    Returns
    -------
    cmap : matplotlib.colors.LinearSegmentedColormap
    """
    data = np.loadtxt(cmap_data_dir / fname, delimiter=',')
    if data.shape[1] != 3:
        raise RuntimeError(f'RGB data files must have 3 columns (got {data.shape[1]})')
    cdict = create_cdict(data[:, 0], data[:, 1], data[:, 2])
    return colors.LinearSegmentedColormap(name, cdict)


def create_cdict(r, g, b):
    """
    Create the color tuples in the correct format.
    """
    i = np.linspace(0, 1, r0.size)

    cdict = {name: list(zip(i, el / 255.0, el / 255.0))
             for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
    return cdict


@u.quantity_input
def suvi_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SUVI images.
    SUVI uses AIA color tables.
    """
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

    # Now create the color dictionary in the correct format
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(
        'GOES-R SUVI {:s}'.format(str(wavelength)), cdict)
