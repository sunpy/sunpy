import re
import json
import warnings

import numpy as np
from packaging.version import Version

import astropy
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import Column, MaskedColumn
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.masked import Masked

import sunpy.coordinates
from sunpy import log
from sunpy.extern import parse
from sunpy.time import parse_time

__all__ = [
    '_freeze',
    '_map_columns_to_times',
    '_map_columns_to_quantities',
    '_map_event_coord_columns_to_coordinates',
    '_map_chain_code_columns_to_coordinates',
]


def _freeze(obj):
    """
    Create a hashable representation of a dict or list.
    Otherwise returns the original input.
    """
    if isinstance(obj, dict):
        return tuple((k, _freeze(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return tuple(_freeze(elem) for elem in obj)
    return obj


def _get_unit_attributes():
    """
    Returns the attributes from the HEK that are used to parse the values in the HEK table.
    """
    unit_properties_filename = get_pkg_data_filename('data/unit_properties.json',
                                                     package='sunpy.net.hek')
    with open(unit_properties_filename) as unit_file:
        unit_properties = json.load(unit_file)
    return unit_properties["attributes"]


def _get_coord_attributes():
    coord_properties_filename = get_pkg_data_filename('data/coord_properties.json',
                                                      package='sunpy.net.hek')
    with open(coord_properties_filename) as coord_file:
        coord_properties = json.load(coord_file)
    return coord_properties["attributes"]


def _map_columns_to_times(table):
    """
    Represent time columns in an HEK response as astropy Time objects.
    """
    # All time columns from https://www.lmsal.com/hek/VOEvent_Spec.html
    time_keys = ['event_endtime', 'event_starttime', 'event_peaktime']
    for tkey in time_keys:
        if tkey in table.colnames:
            table[tkey] = [
                (parse_time(time, format='iso') if time else np.ma.masked) for time in table[tkey]
            ]


def _map_columns_to_quantities(table):
    """
    For columns in an HEK response which have associated units, combine
    data and unit columns into masked columns with units.
    """
    attributes = _get_unit_attributes() + _get_coord_attributes()
    missing_values = {
        'integer': -999,
        'long': -999,
        'float': np.nan,
    }
    dtype_aliases = {
        'string': 'str',
        'integer': 'int',
    }
    for attr in attributes:
        if (name := attr["name"]) not in table.colnames:
            continue
        if attr.get("is_unit_prop", False):
            continue
        if attr.get("is_chaincode", False):
            continue
        mask = np.array([r is None for r in table[name]])
        # NOTE: Fill with missing values explicitly because None cannot be cast to all dtypes.
        data = np.where(mask, missing_values.get(attr['type'], None), table[name])
        if unit_prop := attr.get("unit_prop", False):
            if unit_prop not in table.colnames:
                log.debug(f"Missing unit property {unit_prop} for {name}. Using event_coordunit.")
                unit_prop = "event_coordunit"
            units = np.array([_parse_unit(u_) if u_ is not None else '' for u_ in table[unit_prop]])
            default_unit = u.dimensionless_unscaled if mask.all() else units[~mask][0]
            units = np.where(mask, default_unit, units).tolist()
            # NOTE: This is done per entry because each entry could, in principle,
            # have a different (though compatible) unit.
            data = u.Quantity([d*_u for d, _u in zip(data, units)])
        dtype = dtype_aliases.get(attr['type'], attr['type'])
        table[name] = MaskedColumn(data=data, mask=mask, name=name, dtype=dtype)
    table.remove_columns(
        [attr["name"] for attr in attributes
         if attr.get("is_unit_prop", False) and attr["name"] in table.colnames]
    )


def _map_event_coord_columns_to_coordinates(table):
    """
    For columns in an HEK response which represent the event coordinates, combine
    these columns into a single column with a SkyCoord object.

    These columns are a special case because the event coordinate information is spread across
    many columns and as such multiple columns have to be used to create a single new
    coordinate column.
    """
    # TODO: Confirm that this is the right thing to do. Previous iterations of this
    # code set the frame to ICRS if the event_coordunit was anything but arcsec.
    # See https://www.lmsal.com/hek/VOEvent_Spec.html for more information.
    frame_mapping = {
        "UTC-HGS-TOPO": sunpy.coordinates.HeliographicStonyhurst,
        "UTC-HPC-TOPO": sunpy.coordinates.Helioprojective,
        "UTC-HGC-TOPO": sunpy.coordinates.HeliographicCarrington,
        "UTC-HCR-TOPO": sunpy.coordinates.Heliocentric,
        "UTC-HRC-TOPO": sunpy.coordinates.Heliocentric,  # Possibly a misspelling of HCR?
    }
    event_coords = []
    for row in table:
        data = [row["event_coord1"], row["event_coord2"]]
        if row["event_coord3"] is not None:
            data.append(row["event_coord3"])
        # NOTE: "event_coordunit" can be space or comma-separated string representing the different
        # units of the different coordinate columns or just a single unit string.
        coord_unit = [_parse_unit(unit_string) for unit_string in re.split(r',\s*|\s+', row["event_coordunit"])]
        if len(coord_unit) == 1:
            coord_unit = len(data) * coord_unit
        data = [d*_u for d, _u in zip(data, coord_unit)]
        obstime = row["event_starttime"]
        frame_type = frame_mapping[row["event_coordsys"]]
        frame_kwargs = {"obstime": obstime}
        if frame_type.name != "heliographic_stonyhurst":
            frame_kwargs["observer"] = sunpy.coordinates.get_earth(obstime)
        frame = frame_type(**frame_kwargs)
        representation_type = "spherical"
        # TODO: Confirm whether this is the right way to represent heliocentric radial coordinates or whether
        # this is even what "heliocentric" refers to.
        if isinstance(frame, sunpy.coordinates.Heliocentric):
            representation_type = "cylindrical"
            # NOTE: See entry for Heliocentric Radial in this table:
            # https://docs.sunpy.org/en/stable/reference/coordinates/index.html#supported-coordinate-systems
            data[0] += 90*u.deg
            # The HCC frame expects the data in the reverse order that the HEK returns them.
            data = data[::-1]
            # NOTE: There seem to be cases where event_coord3 is missing for the case of a Heliocentric frame
            if len(data) == 2:
                data.append(1*u.R_sun)
        event_coord = SkyCoord(*data, frame=frame, representation_type=representation_type)
        event_coords.append(event_coord)
    table.add_column(Column(data=event_coords, name="event_coord"))
    # NOTE: Explicitly not removing event_coordunit because it is used as a unit when parsing
    # the other coordinate columns.
    table.remove_columns([
        'event_coordsys',
        'event_coord1',
        'event_coord2',
        'event_coord3',
    ])


def _map_chain_code_columns_to_coordinates(table):
    """
    Convert columns which contain an HEK "chain code", representing either
    a point or a region, into SkyCoord objects.

    In the case of a point, the whole column can be represented as a single SkyCoord
    object with a time varying frame. In the case of a region, the column can be
    represented as a single SkyCoord if each region has the same shape. This is true
    for the bounding box ("bbox") columns. If the shapes are not equal, the column has
    to be expressed as a list of SkyCoord objects.
    """
    # FIXME: There has to be a better way to do this
    frame_class_mapping = {
        'helioprojective': sunpy.coordinates.Helioprojective,
        'heliographic_stonyhurst': sunpy.coordinates.HeliographicStonyhurst,
        'heliocentric': sunpy.coordinates.Heliocentric,
        'heliographic_carrington': sunpy.coordinates.HeliographicCarrington,
        'icrs': ICRS,
    }
    frame_unit_mapping = {
        'helioprojective': [u.arcsec, u.arcsec],
        'heliographic_stonyhurst': [u.deg, u.deg],
        'heliocentric': [u.R_sun, u.deg],
        'heliographic_carrington': [u.deg, u.deg],
    }
    for attr in _get_coord_attributes():
        if not attr.get('is_chaincode', False):
            continue
        if (name := attr["name"]) not in table.colnames:
            continue
        is_point = attr.get('is_point', False)
        frame_kwargs = {}
        if attr['frame'] != 'icrs':
            frame_kwargs['obstime'] = table['event_starttime']
        if attr['frame'] not in ['heliographic_stonyhurst', 'icrs']:
            frame_kwargs['observer'] = 'earth'
        if attr['frame'] == 'heliocentric':
            frame_kwargs['representation_type'] = 'cylindrical'
        frame = frame_class_mapping[attr['frame']](**frame_kwargs)
        coord_data = []
        shape = (1, 2)  # Only used if all entries are masked
        for row in table[name]:
            if row == '' or row is None:
                data = None
            elif is_point:
                data = np.array(parse('POINT({})', row)[0].split(), dtype=float)
                shape = data.shape
            else:
                data = parse('POLYGON(({}))', row)[0]
                data = np.array([r.split() for r in data.split(',')], dtype=float)
                shape = data.shape
            coord_data.append(data)
        # NOTE: Units are cast to arrays of shape (n_rows,2) to make broadcasting
        # to data arrays easier.
        if attr['frame'] in frame_unit_mapping:
            units = np.array([frame_unit_mapping[attr['frame']]]*len(coord_data))
        else:
            units = np.array(list(map(_parse_unit, table['event_coordunit'])))
            units = np.repeat(units[:, np.newaxis], 2, axis=1)
        # NOTE: Filling in masked values with data of appropriate shape allows for
        # broadcasting of coordinate frame information later on if all shapes are the
        # same.
        coord_data = [np.full(shape, np.nan) if c is None else c for c in coord_data]
        if all([coord_data[0].shape==c.shape for c in coord_data]):
            # NOTE: Taking the transpose of the coordinate such that the first
            # dimension of the coordinate is the number of rows in the table. The
            # transpose is not taken of the data prior to construction so that the
            # frames are appropriately broadcast to the data.
            coordinates = _build_masked_coordinate(coord_data, frame, units).T
        else:
            # NOTE: If the shape of the data in each row is not equal, the column has
            # to be expressed as a list of SkyCoord objects and will have a dtype of
            # object.
            if not frame.shape:
                frame = len(coord_data) * [frame]
            coordinates = [_build_masked_coordinate(c, frame[i], units[i]) for i,c in enumerate(coord_data)]
            # FIXME: The mask is created this way for astropy versions below v7 because SkyCoords cannot
            # do not have a mask attribute until v7 and above. Once our minimum version of astropy is v7
            # this conditional can be removed
            if Version(astropy.__version__) < Version('7'):
                mask = [np.isnan(c.cartesian.xyz).all() for c in coordinates]
            else:
                mask = [c.mask.all() for c in coordinates]
            # NOTE: This is expressed explicitly as a column to avoid stacking the coordinates along the
            # dimension that corresponds to the chaincode.
            coordinates = MaskedColumn(data=coordinates,
                                       name=name,
                                       mask=mask,
                                       dtype='object')
        table[name] = coordinates


def _build_masked_coordinate(data, frame, unit):
    # NOTE: Take the transpose so that entries in the list correspond to
    # components of the coordinate.
    coord_data = []
    for _d, _u in zip(np.array(data).T, unit.T):
        # NOTE: This complexity is to allow for broadcasting of units in cases where
        # there is a single unit or an array of units against a data array that may
        # be multidimensional.
        data = np.stack((_d * _u).ravel()).reshape(_d.shape)
        # FIXME: This conditional is because SkyCoords with masked data are only
        # support in astropy v7 and above. Once our minimum version of astropy
        # is v7, this can be removed.
        if Version(astropy.__version__) < Version('7'):
            coord_data.append(data)
        else:
            coord_data.append(Masked(data, mask=np.isnan(_d)))
    if frame.name == 'heliocentric':
        # There is a 90 degree offset between Heliocentric radial and the
        # cylindrical representation of the Heliocentric as defined in sunpy.
        # See https://docs.sunpy.org/en/stable/reference/coordinates/index.html#supported-coordinate-systems.
        coord_data[-1] += 90*u.deg
        # The HEK only returns two coordinates, but a Heliocentric frame requires three.
        # Here, we assume that z=1 R_sun.
        coord_data += [1*u.R_sun]
    return SkyCoord(*coord_data, frame=frame)


def _parse_unit(unit_string):
    """
    Parses an HEK string representation of a unit and converts it into a `astropy.units.Unit`.

    Parameters
    ----------
    unit_string : `str`
        String representation of unit as return by the HEK.

    Returns
    -------
    unit : `~astropy.units.Unit`

    Notes
    ----
    `A complete list of HEK parameters.
    <https://www.lmsal.com/hek/VOEvent_Spec.html>`__
    """
    # NOTE: This has to be done as a regular expression because "ma" is already
    # an astropy unit representing "milliannum" (1/1000 years) and as such it
    # cannot be used as an alias for "milliampere".
    unit_string = re.sub(r"\bma\b", "mA", unit_string)
    # NOTE: These are done separately because unit names with spaces cannot simply be aliased
    # and thus have to be defined as new units. Units without spaces that have unconventional
    # spellings can simply be aliased.
    aliases = {
        "arcseconds": u.arcsec,
        "degrees": u.deg,
        "sec": u.s,
        "emx": u.Mx,
        "amperes": u.A,
        "ergs": u.erg,
        "radians": u.rad,
        "hz": u.Hertz,
        "rsun": u.R_sun,
    }
    enabled_units = [
        u.def_unit("cubic centimeter", represents=u.cm**3),
        u.def_unit("square centimeter", represents=u.cm**2),
        u.def_unit("cubic meter", represents=u.m**3),
        u.def_unit("square meter", represents=u.m**2),
        u.def_unit("ergs per cubic centimeter", represents=u.erg/u.cm**3),
        u.def_unit("ergs/cm^3", represents=u.erg/u.cm**3),
        u.def_unit("HMI pixels", represents=0.5*u.arcsec/u.pixel),
        u.def_unit("HMI pixel", represents=0.5*u.arcsec/u.pixel),
    ]
    # NOTE: This logic accounts for the scenario where multiple units may be contained in
    # a single string. However, this logic should not be applied to custom units which contain
    # spaces, e.g. those custom units enabled above.
    if unit_string not in aliases and unit_string not in [eu.name for eu in enabled_units]:
        split_unit_string = re.split(r',\s*|\s+', unit_string)
        if len(split_unit_string) > 1:
            log.debug(f"Multiple units found in {unit_string}. Using only the first unit.")
            unit_string = split_unit_string[0]
    with u.add_enabled_units(enabled_units), u.set_enabled_aliases(aliases), warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=u.UnitsWarning, message='.*contains multiple slashes*')
        parsed_unit = u.Unit(unit_string)
        return parsed_unit
