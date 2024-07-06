import os
import re
import json
import warnings
from pathlib import Path

import numpy.ma
from regions import PointSkyRegion, PolygonSkyRegion

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Column

from sunpy.time import parse_time

__all__ = [
    'freeze',
    'parse_times',
    'parse_values_to_quantities',
    'UNIT_FILE_PATH',
    'COORD_FILE_PATH'
]

UNIT_FILE_PATH = Path(os.path.dirname(__file__)) / "unit_properties.json"
COORD_FILE_PATH = Path(os.path.dirname(__file__)) / "coord_properties.json"

def freeze(obj):
    """
    Create a hashable representation of a dict or list.
    Otherwise returns the original input.
    """
    if isinstance(obj, dict):
        return tuple((k, freeze(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return tuple(freeze(elem) for elem in obj)
    return obj

# NOTE: Needs unit test
def parse_times(table):
    """
    Parses the time values in an Astropy table into `astropy.time.Time`.

    Parameters
    ----------
    table : `astropy.table.Table`
        Astropy table.

    Returns
    -------
    `astropy.table.Table`
    """
    # All time columns from https://www.lmsal.com/hek/VOEvent_Spec.html
    time_keys = ['event_endtime', 'event_starttime', 'event_peaktime']
    for tkey in time_keys:
        if tkey in table.colnames:
            table[tkey] = [
                (parse_time(time, format='iso') if time else numpy.ma.masked) for time in table[tkey]
            ]
    return table

def parse_values_to_quantities(table):
    """
    Parses the values in an Astropy table into Astropy objects.

    Parameters
    ----------
    table : `astropy.table.Table`
        Astropy table.

    Returns
    -------
    `astropy.table.Table`

    Raises
    ------
    TypeError
        If `table` is not an Astropy table.
    """
    with open(UNIT_FILE_PATH) as unit_file:
        unit_properties = json.load(unit_file)
    unit_attributes = unit_properties["attributes"]

    with open(COORD_FILE_PATH) as coord_file:
        coord_properties = json.load(coord_file)
    coord_attributes = coord_properties["attributes"]
    table = parse_columns_to_table(table, unit_attributes)
    table = parse_columns_to_table(table, coord_attributes, is_coord_prop = True)
    return table

def parse_columns_to_table(table, attributes, is_coord_prop = False):
    """
    Parses the columns in an Astropy table and convert the values into Astropy objects.

    Parameters
    ----------
    table: `astropy.table.Table`
        Astropy table.
    attributes : list
        A list of HEK unit attributes or coordinate attributes.
    is_coord_prop : bool
        To specify if ``attributes`` is a list of unit attributes or coordinate attributes.

    Returns
    -------
    `astropy.table.Table`

    Raises
    ------
    TypeError
        If ``table`` is not an Astropy table.
    KeyError
        If any of the attribute dictionaries are missing required keys (i.e. "name", "unit_prop").
    """
    if is_coord_prop:
        event_coord_col = []
        for idx in range(len(table['event_coord1'])):
            unit = get_unit(table['event_coordunit'][idx])
            coord1 = table['event_coord1'][idx]
            coord2 = table['event_coord2'][idx]
            coord3 = table['event_coord3'][idx]
            frame = 'helioprojective' if unit==u.arcsec else 'icrs'
            if coord3:
                event_coord = SkyCoord(coord1*unit, coord2*unit, coord3*unit, frame=frame)
            else:
                event_coord = SkyCoord(coord1*unit, coord2*unit, frame=frame)
            event_coord_col.append(event_coord)
        event_coord_col = Column(event_coord_col, name='event_coord')
        del table['event_coord1']
        del table['event_coord2']
        del table['event_coord3']
        table.add_column(event_coord_col)

    for attribute in attributes:
        if attribute["name"] in table.colnames and ("unit_prop" in attribute or attribute.get("is_chaincode", False)) and attribute.get("is_unit_prop", True):
            unit_attr = ""
            if is_coord_prop:
                unit_attr = "event_coordunit"
            else:
                unit_attr = attribute["unit_prop"]

            new_column = []
            for idx, value in enumerate(table[attribute["name"]]):
                new_value = ""
                if value in ["", None]:
                    new_value = value
                elif attribute.get("is_chaincode", False):
                    unit = table[attribute["unit_prop"]][idx]
                    time = table['event_starttime'][idx]
                    new_value = parse_chaincode(value, attribute, unit, time)
                else:
                    unit = get_unit(table[unit_attr][idx])
                    new_value = value * unit
                    if attribute["name"] == 'obs_meanwavel':
                        new_value = new_value.to(u.angstrom)
                new_column.append(new_value)

            if not attribute.get("is_chaincode"):
                new_column = Column(new_column, name=table[attribute["name"]], dtype=u.Quantity)
            table[attribute["name"]] = new_column

    for attribute in attributes:
        if attribute.get("is_unit_prop", False) and attribute["name"] in table.colnames:
            del table[attribute["name"]]
    return table

def parse_chaincode(value, attribute, unit, time):
    """
    Parses a string representation of coordinates and convert them into a PolygonSkyRegion object
    using units based on the specified coordinate frame.

    Parameters
    ----------
    value : str
        A polygon defined using vertices in sky coordinates.
    attribute : dict
        An object from "coord_properties.json"
    unit : str
        The unit of the coordinates
    time : `~astropy.time.core.Time`
        An event_starttime row parsed into astropy time.

    Returns
    -------
    `PolygonSkyRegion`
        A polygon defined using vertices in sky coordinates.

    Raises
    ------
    IndexError
        Because ``value`` does not contain the expected '((' and '))' substrings.
    UnitConversionError
        Because the units set by ``coord1_unit`` or ``coord2_unit`` are incompatible with the values being assigned.
    """
    observer = 'earth' # There is an assumption that earth is the observer.
    coord1_unit = u.deg
    coord2_unit = u.deg
    if attribute["frame"] == "helioprojective":
        coord1_unit = u.arcsec
        coord2_unit = u.arcsec
    elif attribute["frame"] == "heliocentric":
        coord1_unit = u.R_sun # Nominal solar radius
    elif attribute["frame"] == "icrs":
        coord1_unit = get_unit(unit)
        coord2_unit = get_unit(unit)

    if attribute.get("is_point"):
        coordinates = value.replace("POINT(", "").replace(")", "").split()
        coord_list = [float(coordinate) for coordinate in coordinates]
        coord_list[0] *= coord1_unit
        coord_list[1] *= coord2_unit
        if attribute["frame"] == "heliocentric":
            center_sky = SkyCoord(coord_list[0], coord_list[1], [1]* len(coord_list) * u.AU, obstime=time, observer=observer, representation_type="cylindrical", frame=attribute["frame"])
            return PolygonSkyRegion(vertices=center_sky)
        else:
            center_sky = SkyCoord(coord_list[0], coord_list[1], obstime=time, observer=observer, frame=attribute["frame"])
        return PointSkyRegion(center=center_sky)
    coordinates_str = value.split('((')[1].split('))')[0]
    coord1_list = [float(coord.split()[0]) for coord in coordinates_str.split(',')] * coord1_unit
    coord2_list = [float(coord.split()[1]) for coord in coordinates_str.split(',')] * coord2_unit
    vertices = {}
    if attribute["frame"] == "heliocentric":
        vertices = SkyCoord(coord1_list, coord2_list, [1]* len(coord1_list) * u.AU, obstime=time, observer=observer, representation_type="cylindrical", frame="heliocentric")
    else:
        vertices = SkyCoord(coord1_list, coord2_list, obstime=time, observer=observer, frame=attribute["frame"])

    return PolygonSkyRegion(vertices=vertices)

def get_unit(unit):
    """
    Converts string into astropy unit.

    Parameters
    ----------
    unit : str
        The targeted unit

    Returns
    -------
    unit : `~astropy.units.Unit`
        Astropy unit object (e.g. :`~astropy.units.Unit` or :`~astropy.units.CompositeUnit`)

    Raises
    ------
    ValueError
        Because ``unit`` did not parse as unit.

    Notes
    ----
    `A complete list of HEK parameters.
    <https://www.lmsal.com/hek/VOEvent_Spec.html>`__
    """
    cm2 = u.def_unit("cm2", u.cm**3)
    m2 = u.def_unit("m2", u.m**2)
    m3 = u.def_unit("m3", u.m**3)
    erg_per_cm3 = u.def_unit("ergs/cm^3", u.erg/u.ml)

    aliases = {
        "steradian": u.sr,
        "arcseconds": u.arcsec,
        "degrees": u.deg,
        "sec": u.s,
        "emx": u.Mx,
        "amperes": u.A,
        "ergs": u.erg,
        "cubic centimeter": u.ml,
        "square centimeter": cm2,
        "cubic meter": m3,
        "square meter": m2,
        "ergs per cubic centimeter": erg_per_cm3,
    }
    with u.add_enabled_units([cm2, m2, m3]), u.set_enabled_aliases(aliases), warnings.catch_warnings():
        # Units for coordinate frames have more than one unit, otherwise it will be just one unit.
        # There is an assumption that coord1_unit, coord2_unit and coord3_unit are the same.
        warnings.filterwarnings("ignore", category=u.UnitsWarning)
        if unit in aliases:
            unit = u.Unit(aliases[unit])
        else:
            unit = u.Unit(re.split(r'[, ]', unit)[0].lower())
        return unit
