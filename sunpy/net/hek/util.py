import os
import re
import json
from pathlib import Path

from regions import PolygonSkyRegion

from astropy.coordinates import SkyCoord
from astropy import units as u

from sunpy.time import parse_time

UNIT_FILE_PATH = Path(os.path.dirname(__file__)) / "unit_properties.json"
COORD_FILE_PATH = Path(os.path.dirname(__file__)) / "coord_properties.json"

def freeze(obj):
    """ Create hashable representation of result dict. """
    if isinstance(obj, dict):
        return tuple((k, _freeze(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return tuple(_freeze(elem) for elem in obj)
    return obj

def parse_times(table):
    # All time columns from https://www.lmsal.com/hek/VOEvent_Spec.html
    time_keys = ['event_endtime', 'event_starttime', 'event_peaktime']
    for tkey in time_keys:
        if tkey in table.colnames and not any(time == "" for time in table[tkey]):
            table[tkey] = parse_time(table[tkey])
            table[tkey].format = 'iso'
    return table

def parse_unit(table, attribute, is_coord_prop = False):
    if attribute.get("is_chaincode", False):
        return table
    unit_attr=""
    if is_coord_prop:
        unit_attr = "event_coordunit"
    else:
        unit_attr = attribute["unit_prop"]
    for row in table:
        if unit_attr in table.colnames and row[unit_attr] not in ["", None] and table[attribute["name"]].unit is not None:
            table[attribute["name"]].unit = get_unit(attribute["unit_prop"], row[unit_attr], is_coord_prop= is_coord_prop)
            break
    return table

def parse_astropy_unit(str):
    try:
        unit = u.Unit(str)
    except ValueError:
        try:
            unit = u.Unit(str.lower())
        except ValueError:
            try:
                unit = u.Unit(str.capitalize())
            except ValueError:
                units = str.split(" per ")
                unit = None
                if len(units) >1:
                    unit = parse_astropy_unit(units[0])
                    for idx in range(1, len(units)):
                        unit = unit/parse_astropy_unit(units[idx])
                else:
                    unit = unit_mapping[str]

    return unit

def get_unit(unit_prop, str, is_coord_prop = False):
    if is_coord_prop:
        coord1_unit, coord2_unit, coord3_unit = None, None, None
        coord_units = re.split(r'[, ]', str)
        if len(coord_units) == 1: # deg
           coord1_unit = coord2_unit = parse_astropy_unit(coord_units[0])
        elif len(coord_units) == 2:
            coord1_unit = parse_astropy_unit(coord_units[0])
            coord2_unit = parse_astropy_unit(coord_units[1])
        else:
            coord1_unit = parse_astropy_unit(coord_units[0])
            coord2_unit = parse_astropy_unit(coord_units[1])
            coord3_unit = parse_astropy_unit(coord_units[2])
        return locals()[unit_prop]
    else:
        return parse_astropy_unit(str)

def parse_chaincode(value, idx, attribute, unit_prop):
    coord1_unit = u.deg
    coord2_unit = u.deg
    if attribute["frame"] == "helioprojective":
        coord1_unit = u.arcsec
        coord2_unit = u.arcsec
    elif attribute["frame"] == "heliocentric":
        coord1_unit = u.R_sun
        coord2_unit = u.deg
    elif attribute["frame"] == "icrs":
        coord1_unit = get_unit("coord1_unit", unit_prop, is_coord_prop = True)
        coord2_unit = get_unit("coord2_unit", unit_prop, is_coord_prop = True)

    coordinates_str = value.split('((')[1].split('))')[0]
    coord1_list = [float(coord.split()[0]) for coord in coordinates_str.split(',')] * coord1_unit
    coord2_list = [float(coord.split()[1]) for coord in coordinates_str.split(',')] * coord2_unit
    vertices = {}
    if attribute["frame"] == "heliocentric":
       vertices = SkyCoord(coord1_list, coord2_list, [1]* len(coord1_list)* u.AU, representation_type="cylindrical" , frame="heliocentric" )
    else:
        vertices = SkyCoord(coord1_list, coord2_list, frame=attribute["frame"])

    return PolygonSkyRegion(vertices = vertices)

def parse_columns_to_table(table, attributes, is_coord_prop = False):
    for attribute in attributes:
        if attribute.get("is_unit_prop", False):
            pass
        elif attribute["name"] in table.colnames and ("unit_prop" in attribute or attribute.get("is_chaincode", False)):
            table = parse_unit(table, attribute, is_coord_prop)
            unit_attr = ""
            if is_coord_prop:
                unit_attr = "event_coordunit"
            else:
                unit_attr = attribute["unit_prop"]

            new_column = []
            for idx, value in enumerate(table[attribute["name"]]):
                if value in ["", None]:
                    new_column.append(value)
                elif attribute.get("is_chaincode", False):
                    new_column.append(parse_chaincode(value, idx, attribute, table[attribute["unit_prop"]][idx]))
                else:
                    new_column.append(value * get_unit(attribute["unit_prop"], table[unit_attr][idx], is_coord_prop= is_coord_prop))
            table[attribute["name"]] = new_column
    for attribute in attributes:
        if attribute.get("is_unit_prop", False) and attribute["name"] in table.colnames:
            del table[attribute["name"]]
    return table

def parse_values_to_quantities(table):
    with open(UNIT_FILE_PATH) as unit_file:
        unit_properties = json.load(unit_file)
    unit_attributes = unit_properties["attributes"]

    with open(COORD_FILE_PATH) as coord_file:
        coord_properties = json.load(coord_file)
    coord_attributes = coord_properties["attributes"]
    table = parse_columns_to_table(table, unit_attributes)
    table = parse_columns_to_table(table, coord_attributes, is_coord_prop= True)
    return table


