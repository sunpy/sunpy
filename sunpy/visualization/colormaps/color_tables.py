"""
This module provides enums for generating
`~matplotlib.colors.LinearSegmentedColormaps`, and a dictionary of these
dictionaries.
"""
import pathlib

import matplotlib.colors as colors
import numpy as np

import astropy.units as u
import enum

__all__ = [
    '_cmap_from_rgb', 'cmap_from_rgb_file', 'iris_sji_color_table', 'xrt_color_table'
]

cmap_data_dir = pathlib.Path(__file__).parent.absolute() / 'data'

sxt_gold_r = np.concatenate((np.linspace(0, 255, num=185,
                                     endpoint=False), 255 * np.ones(71)))
sxt_gold_g = 255 * (np.arange(256)**1.25) / (255.0**1.25)
sxt_gold_b = np.concatenate((np.zeros(185), 255.0 * np.arange(71) / 71.0))
grayscale = np.arange(256)

idl_3 = np.loadtxt(cmap_data_dir / 'idl_3.csv', delimiter=',')
r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]

c0 = np.arange(256, dtype='f')
c1 = (np.sqrt(c0) * np.sqrt(255.0)).astype('f')
c2 = (np.arange(256)**2 / 255.0).astype('f')
c3 = ((c1 + c2 / 2.0) * 255.0 / (c1.max() + c2.max() / 2.0)).astype('f')

def _cmap_from_rgb(color, name):
    r,g, b = color
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
    cmap : matplotlib.colors.LinearSegmentedColormap
    """
    data = np.loadtxt(cmap_data_dir / fname, delimiter=',')
    if data.shape[1] != 3:
        raise RuntimeError(f'RGB data files must have 3 columns (got {data.shape[1]})')
    return _cmap_from_rgb((data[:, 0],data[:, 1],data[:, 2]), name)
    
def xrt_color_table():
    """
    Returns the color table used for all Hinode XRT images.
    """
    return _cmap_from_rgb(r0, g0, b0, 'Hinode XRT')

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
        data = color_table[measurement]
    except KeyError:
        raise ValueError("Invalid IRIS SJI waveband.  Valid values are \n" +
                         str(list(color_table.keys())))

    return _cmap_from_rgb(data, f'IRIS SJI {measurement:s}')

def create_cdict(r, g, b):
    """
    Create the color tuples in the correct format.
    """
    i = np.linspace(0, 1, r0.size)

    cdict = {name: list(zip(i, el / 255.0, el / 255.0))
             for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
    return cdict

class Colors(enum.Enum):
    
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
    color = {171*u.angstrom: 'dark_blue', 195*u.angstrom: 'dark_green',
                 284*u.angstrom: 'yellow', 304*u.angstrom: 'dark_red'
    }
    sxt_colors = {
            'al': (sxt_gold_r, sxt_gold_g, sxt_gold_b),
            'wh': (grayscale, grayscale, grayscale)
        }
    rgb_colors = {
            'intensity': (r0, g0, b0),
        }