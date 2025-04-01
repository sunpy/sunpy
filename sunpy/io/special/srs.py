"""
This module implements a SRS File Reader.
"""
import re
import datetime
from collections import OrderedDict

import numpy as np

import astropy.io.ascii
import astropy.units as u
from astropy.table import Column, MaskedColumn, QTable, vstack

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

    header, section_lines, supplementary_lines = split_lines(file_lines)

    return make_table(header, section_lines, supplementary_lines)


def make_table(header, section_lines, supplementary_lines):
    """
    From the separated section lines and the header, clean up the data and
    convert to a `~astropy.table.QTable`.
    """
    meta_data = get_meta_data(header, supplementary_lines)

    tables = []
    for i, lines in enumerate(section_lines):
        if lines:
            key = list(meta_data['id'].keys())[i]
            t1 = astropy.io.ascii.read(lines)

            # Change column names into titlecase
            column_names = list(t1.columns)
            t1.rename_columns(column_names, new_names=[col.title() for col in column_names])

            if len(t1) == 0:
                col_data_types = {
                    # ID : <class 'str'>
                    'Nmbr': np.dtype('i4'),
                    'Location': np.dtype('U6'),
                    'Lo': np.dtype('i8'),
                    'Area': np.dtype('i8'),
                    'Z': np.dtype('U3'),
                    'Ll': np.dtype('i8'),
                    'Nn': np.dtype('i8'),
                    'Magtype': np.dtype('S4'),
                    'Lat': np.dtype('i8'),
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
    column_mapping = {
        'Nmbr': 'Number',
        'Nn': 'Number of Sunspots',
        'Lo': 'Carrington Longitude',
        'Magtype': 'Mag Type',
        'Ll': 'Longitudinal Extent',
    }

    for old_name, new_name in column_mapping.items():
        out_table.rename_column(old_name, new_name)

    # Define a Solar Hemisphere Unit
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
    return the lines for the header, a list of lines for each section that
    is not 'None', and a list of supplementary lines after the main sections
    if not 'None'.
    """
    section_lines = []
    final_section_lines = []
    for i, line in enumerate(file_lines):
        if re.match(r'^(I\.|IA\.|II\.)', line):
            section_lines.append(i)
        if re.match(r'^(III|COMMENT|EFFECTIVE 2 OCT 2000|PLAIN|This message is for users of the NOAA/SEC Space|NNN)', line, re.IGNORECASE):
            final_section_lines.append(i)

    if final_section_lines and final_section_lines[0] > section_lines[-1]:
        section_lines.append(final_section_lines[0])

    header = file_lines[:section_lines[0]]
    header += [file_lines[s] for s in section_lines]

    # Append comments to the comment lines
    for line in section_lines:
        file_lines[line] = '# ' + file_lines[line]
    t1_lines = file_lines[section_lines[0]:section_lines[1]]
    # Remove the space so table reads it correctly
    t1_lines[1] = re.sub(r'Mag\s*Type', r'Magtype', t1_lines[1], flags=re.IGNORECASE)
    t2_lines = file_lines[section_lines[1]:section_lines[2]]

    # SRS files before 2000-10-02 files may have an empty `COMMENT` column in ``t2_lines``
    if "COMMENT" in t2_lines[1].split():
        expected_pattern_dict = {
            'Nmbr': r'^\d+$',
            'Location': r'^(?:[NESW](?:\d{2})){1,2}$',
            'Lo': r'^\d+$',
        }
        # Try to drop the comment column and return in original format
        t2_lines[1:] = _try_drop_empty_column("COMMENT", t2_lines[1:], expected_pattern_dict)

    if len(section_lines) > 3:
        t3_lines = file_lines[section_lines[2]:section_lines[3]]
        supplementary_lines = file_lines[section_lines[3]:]
    else:
        t3_lines = file_lines[section_lines[2]:]
        supplementary_lines = None

    lines = [t1_lines, t2_lines, t3_lines]
    for i, ll in enumerate(lines):
        if len(ll) > 2 and ll[2].strip().title() == 'None':
            del ll[2]

    return header, lines, supplementary_lines


def get_meta_data(header, supplementary_lines):
    """
    Convert a list of header lines and a list of supplementary lines (if not 'None') into a meta data dict.
    """
    meta_lines = []
    for line in header:
        if line.startswith(':'):
            meta_lines.append(line)

    meta_data = {}
    for m in meta_lines:
        if re.search(r'Corrected\s*Copy', m, re.IGNORECASE):
            meta_data['corrected'] = True
            continue
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

    if supplementary_lines:
        meta_data['supplementary_lines'] = [sl.strip() for sl in supplementary_lines]

    return meta_data


def parse_longitude(value):
    """
    Parse longitude in the form "W10" or "E10".
    """
    lonsign = {'W': 1, 'E': -1}
    if "W" in value or "E" in value:
        return lonsign[value[3]] * float(value[4:])


def parse_latitude(value):
    """
    Parse latitude in the form "S10" or "N10".
    """
    latsign = {'N': 1, 'S': -1}
    if "N" in value or "S" in value:
        return latsign[value[0]] * float(value[1:3])


def parse_location(column):
    """
    Given a column of location data in the form "S10E10" convert to two columns
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
    Given an input column of "latitudes" in the form "S10" parse them and add
    them to an existing column of "latitudes".
    """
    for i, loc in enumerate(column):
        if loc:
            latitude_column.mask[i] = False
            latitude_column[i] = parse_latitude(loc)
    return latitude_column


def _try_drop_empty_column(column_name_to_drop, data_lines, pattern_dict):
    """
    Try dropping an empty ``column_name_to_drop`` from ``data_lines``.

    Parameters
    ----------
    column_name_to_drop : `str`
        Name of the empty column to be dropped.
    data_lines : `list[str]`
        List of lines extracted from a file (each line is a string)
        corresponding to the header (e.g. ``header = data_lines[0]``)
        and the data (``data = data_lines[1:]``)
    pattern_dict : `dict`
        A dictionary specifying the patterns to match for each column

    Returns
    -------
    `list[str]`
        The modified ``data_lines`` in titlecase with the specified column dropped, if all validations pass.

    """
    # Create a lowercase pattern dict
    pattern_dict_lower = {key.lower(): value for key, value in pattern_dict.items()}

    # Extract columns and rows
    header_line, *row_lines = data_lines
    column_list = [column.strip().lower() for column in header_line.split()]

    # Drop ``column_name_to_drop`` if exists
    try:
        column_index = column_list.index(column_name_to_drop.strip().lower())
        column_list.pop(column_index)
    except ValueError:
        raise ValueError(f"The column '{column_name_to_drop}' does not exist.")

    # Remove the dropped column from pattern_dict
    pattern_dict_lower.pop(column_name_to_drop.strip().lower(), None)

    # If the data is `None`, just return the header/data
    if row_lines[0].strip().title() == 'None':
        # Return as titlecase
        column_list = [col.title() for col in column_list]
        return [" ".join(column_list)] + row_lines

    # Check if the remaining columns are a subset of the columns in pattern_dict
    remaining_columns_set = set(column_list)
    pattern_columns_set = set(pattern_dict_lower.keys())
    if not remaining_columns_set.issubset(pattern_columns_set):
        raise ValueError("The remaining columns are not a subset of the columns in ``pattern_dict``.")

    # Check if all rows have the same length as the remaining columns
    row_lengths_equal = all(len(row.split()) == len(column_list) for row in row_lines)
    if not row_lengths_equal:
        raise ValueError("not all rows have the same number of values as the remaining columns.")

    # Check that the row values are consistent with the provided pattern dictionary
    matching_pattern = all(all(re.match(pattern_dict_lower[column], value) for column, value in zip(column_list, row.split())) for row in row_lines)
    if not matching_pattern:
        raise ValueError("not all rows match the provided pattern.")

    # Return as titlecase
    column_list = [col.title() for col in column_list]
    return [" ".join(column_list)] + row_lines
