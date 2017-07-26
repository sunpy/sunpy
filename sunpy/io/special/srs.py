"""
This module implements SRS File Reader.
"""
import datetime
from collections import OrderedDict

import numpy as np
from astropy.table import QTable, MaskedColumn, Column, vstack
import astropy.io.ascii
import astropy.units as u

__all__ = ['read_srs']


def read_srs(filepath):
    """
    Parse a SRS table from NOAA SWPC.

    Parameters
    ----------
    filepath : `str`
        The full path to a SRS table.

    Returns
    -------

    table : `astropy.table.QTable`
        Table containing a stacked table from all the tables in the SRS file.
        The header information is stored in the ``.meta`` attribute.

    """
    with open(filepath) as srs:
        file_lines = srs.readlines()

    header, section_lines = split_lines(file_lines)

    return make_table(header, section_lines)


def make_table(header, section_lines):
    """
    From the seperated section lines and the header, clean up the data and
    convert to a QTable.
    """
    meta_data = get_meta_data(header)

    tables = []
    for i, lines in enumerate(section_lines):
        if lines:
            key = list(meta_data['id'].keys())[i]
            t1 = astropy.io.ascii.read(lines)

            if len(t1) == 0:
                col_data_types = {
                    # ID : <class 'str'>
                    'Nmbr': np.dtype('i4'),
                    'Location': np.dtype('U6'),
                    'Lo': np.dtype('i8'),
                    'Area': np.dtype('i8'),
                    'Z': np.dtype('U3'),
                    'LL': np.dtype('i8'),
                    'NN': np.dtype('i8'),
                    'MagType': np.dtype('S4'),
                    'Lat': np.dtype('i8')
                }
                for c in t1.itercols():
                    # Put data types of columns in empty table to correct types,
                    # or else vstack will fail.
                    c.dtype = col_data_types[c._name]
                t1.add_column(
                    Column(data=None, name="ID", dtype=('S2')), index=0)
            else:
                t1.add_column(Column(data=[key] * len(t1), name="ID"), index=0)

            tables.append(t1)

    out_table = vstack(tables)

    # Parse the Location column in Table 1
    if 'Location' in out_table.columns:
        col_lat, col_lon = parse_location(out_table['Location'])
        del out_table['Location']
        out_table.add_column(col_lat)
        out_table.add_column(col_lon)

    # Parse the Lat column in Table 3
    if 'Lat' in out_table.columns:
        parse_lat_col(out_table['Lat'], out_table['Latitude'])
        del out_table['Lat']

    # Give columns more sensible names
    out_table.rename_column("Nmbr", "Number")
    out_table.rename_column("NN", "Number of Sunspots")
    out_table.rename_column("Lo", "Carrington Longitude")
    out_table.rename_column("MagType", "Mag Type")
    out_table.rename_column("LL", "Longitudinal Extent")

    # Define a Solar Hemispere Unit
    a = {}
    u.def_unit(
        "SH",
        represents=(2 * np.pi * u.solRad**2),
        prefixes=True,
        namespace=a,
        doc="A solar hemisphere is the area of the visible solar disk.")

    # Set units on the table
    out_table['Carrington Longitude'].unit = u.deg
    out_table['Area'].unit = a['uSH']
    out_table['Longitudinal Extent'].unit = u.deg

    out_table.meta = meta_data

    # Number should be formatted in 10000 after 2002-06-15.
    if out_table.meta['issued'] > datetime.datetime(2002, 6, 15):
        out_table['Number'] += 10000

    return QTable(out_table)


def split_lines(file_lines):
    """
    Given all the lines in the file split based on the three sections and
    return the lines for the header and a list of lines for each section that
    is not 'None'
    """
    section_lines = []
    for i, line in enumerate(file_lines):
        if line.startswith(("I.", "IA.", "II.")):
            section_lines.append(i)

    header = file_lines[:section_lines[0]]
    header += [file_lines[s] for s in section_lines]

    # Append comments to the comment lines
    for l in section_lines:
        file_lines[l] = '# ' + file_lines[l]

    t1_lines = file_lines[section_lines[0]:section_lines[1]]
    # Remove the space so table reads it correctly
    t1_lines[1] = t1_lines[1].replace('Mag Type', 'MagType')
    t2_lines = file_lines[section_lines[1]:section_lines[2]]
    t3_lines = file_lines[section_lines[2]:]

    lines = [t1_lines, t2_lines, t3_lines]
    for i, ll in enumerate(lines):
        if ll[2].strip() == 'None':
            del ll[2]

    return header, lines


def get_meta_data(header):
    """
    Convert a list of header lines into a meta data dict.
    """
    meta_lines = []
    for l in header:
        if l.startswith(':'):
            meta_lines.append(l)

    meta_data = {}
    for m in meta_lines:
        k, v = m.strip().split(':')[1:]
        meta_data[k.lower()] = v.strip()
    meta_data['issued'] = datetime.datetime.strptime(meta_data['issued'],
                                                     "%Y %b %d %H%M UTC")

    # Get ID descriptions
    meta_data['id'] = OrderedDict()
    for h in header:
        if h.startswith(("I.", "IA.", "II.")):
            i = h.find('.')
            k = h[:i]
            v = h[i + 2:]
            meta_data['id'][k] = v.strip()

    meta_data['header'] = [h.strip() for h in header]
    return meta_data


def parse_longitude(value):
    """
    Parse longitude in the form 'W10' or 'E10'
    """
    lonsign = {'W': 1, 'E': -1}
    if "W" in value or "E" in value:
        return lonsign[value[3]] * float(value[4:])


def parse_latitude(value):
    """
    Parse latitude in the form 'S10' or 'N10'
    """
    latsign = {'N': 1, 'S': -1}
    if "N" in value or "S" in value:
        return latsign[value[0]] * float(value[1:3])


def parse_location(column):
    """
    Given a column of location data in the form 'S10E10' convert to two columns
    of angles.
    """
    latitude = MaskedColumn(name="Latitude", unit=u.deg)
    longitude = MaskedColumn(name="Longitude", unit=u.deg)

    for i, loc in enumerate(column):
        if loc:
            lati = parse_latitude(loc)
            longi = parse_longitude(loc)
            latitude = latitude.insert(i, lati)
            longitude = longitude.insert(i, longi)
        else:
            latitude = latitude.insert(i, None, mask=True)
            longitude = longitude.insert(i, None, mask=True)
    return latitude, longitude


def parse_lat_col(column, latitude_column):
    """
    Given an input column of Latitudes in the form 'S10' parse them and add
    them to an existing column of Latitudes.
    """
    for i, loc in enumerate(column):
        if loc:
            latitude_column.mask[i] = False
            latitude_column[i] = parse_latitude(loc)
    return latitude_column
