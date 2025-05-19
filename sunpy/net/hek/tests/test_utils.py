import copy

import numpy as np
import pytest

import astropy.table
import astropy.time
import astropy.units as u
import astropy.utils
import astropy.utils.masked
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from sunpy.coordinates import (
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    get_earth,
)
from sunpy.net.hek.utils import (
    _map_chain_code_columns_to_coordinates,
    _map_columns_to_quantities,
    _map_event_coord_columns_to_coordinates,
    _parse_unit,
)


@pytest.mark.parametrize(('input_unit', 'expected_unit'), [
    ('DN/sec/pixel', u.DN/(u.pix*u.s)),
    ('ergs per cubic centimeter', u.erg/u.cm**3),
    ('HMI pixel', 0.5*u.arcsec/u.pixel),
    ('HMI pixels', 0.5*u.arcsec/u.pixel),
    ('ma/m2', u.milliampere/u.m**2),
    ('deg deg', u.deg),
    ('arcsec, arcsec', u.arcsec),
    ('arcsec,arcsec', u.arcsec),
    ("arcseconds", u.arcsec),
    ("degrees", u.deg),
    ("sec", u.s),
    ("emx", u.Mx),
    ("amperes", u.A),
    ("ergs", u.erg),
    ("radians", u.rad),
    ("hz", u.Hertz),
    ("rsun", u.R_sun),
    ("cubic centimeter", u.cm**3),
    ("square centimeter", u.cm**2),
    ("cubic meter", u.m**3),
    ("square meter", u.m**2),
    ("ergs/cm^3", u.erg/u.cm**3),
])
def test_parse_unit(input_unit, expected_unit):
    assert _parse_unit(input_unit) == expected_unit


def test_parse_unit_invalid():
    """
    Test that an unexpected unit errors in the same way as u.Unit
    """
with pytest.raises(ValueError, match=r"'foo' did not parse as unit"):
    _parse_unit('foo')


def map_event_coords_to_column(coord):
    """
    Given a coordinate, create a table that mimics what is returned by the HEK. This is to test
    roundtripping of coordinates through the logic used to parse the "event_coord" columns.
    """
    frame_sys_mapping = {
        "heliographic_stonyhurst": "UTC-HGS-TOPO",
        "helioprojective": "UTC-HPC-TOPO",
        "heliographic_carrington": "UTC-HGC-TOPO",
        "heliocentric": "UTC-HRC-TOPO",
    }
    columns = {
        "event_coordsys": [frame_sys_mapping[coord.frame.name]],
        "event_coordunit": [coord.info.unit],
        "event_starttime": [coord.obstime.iso],
    }
    for i, c in enumerate(coord.data.components):
        columns[f"event_coord{i+1}"] = [getattr(coord.data, c).value]
    # Swap order of first two coords and delete third in the case of HCC (radial)
    if coord.frame.name == "heliocentric":
        columns["event_coord1"], columns["event_coord2"] = columns["event_coord2"], columns["event_coord1"]
        columns["event_coord1"] = [columns["event_coord1"][0] - 90]
        del columns["event_coord3"]
        columns['event_coordunit'] = [', '.join(columns["event_coordunit"][0].split(",")[:-1][::-1])]
    return astropy.table.Table(columns)


def test_event_coord_column_mapping():
    obstime = "2020-01-01"
    observer = get_earth(obstime)
    input_coordinates = [
        SkyCoord(100,200, unit='arcsec', frame=Helioprojective(observer=observer)),
        SkyCoord(20, 45, unit='degree', frame=HeliographicStonyhurst(obstime=obstime)),
        SkyCoord(20*u.deg, 45*u.deg, 1*u.R_sun, frame=HeliographicStonyhurst(obstime=obstime)),
        SkyCoord(180*u.deg, 0*u.deg, frame=HeliographicCarrington(observer=observer)),
        SkyCoord(0.1, 45, 1, unit=['R_sun', 'deg', 'R_sun'], frame=Heliocentric(observer=observer), representation_type='cylindrical'),
    ]
    hek_table = astropy.table.vstack([map_event_coords_to_column(ic) for ic in input_coordinates])
    _map_event_coord_columns_to_coordinates(hek_table)
    assert "event_coordunit" in hek_table.colnames
    for ic, hc in zip(input_coordinates, hek_table['event_coord']):
        assert ic == hc


def test_quantity_column_mapping():
    data = {
        "meanvertcurrentdensity": [1.0, None, 3.0],
        "currentdensityunit": ["ma/m2", None, "ma/m2"],
        "meanphotoenergydensity": [20.0, None, 30.0],
        "meanenergydensityunit": ["ergs per cubic centimeter", None, "ergs per cubic centimeter"],
        "boundbox_c1ur": [200, 250, 0],
        "event_coordunit": ["arcsec, arcsec", "arcsec", "deg deg"],
    }
    hek_table = astropy.table.Table(data)
    _map_columns_to_quantities(hek_table)
    assert u.allclose(hek_table["meanvertcurrentdensity"],
                      astropy.utils.masked.Masked(data["meanvertcurrentdensity"]*u.mA/u.m**2, mask=[False, True, False]))
    assert u.allclose(hek_table["meanphotoenergydensity"],
                      astropy.utils.masked.Masked(data["meanphotoenergydensity"]*u.erg/u.cm**3, mask=[False, True, False]))
    assert u.allclose(hek_table["boundbox_c1ur"], data["boundbox_c1ur"]*u.arcsec)
    assert hek_table.colnames == ["meanvertcurrentdensity", "meanphotoenergydensity", "boundbox_c1ur"]


obstime = astropy.time.Time("2010-01-01") + np.arange(5) * u.year
hpc_frame = Helioprojective(observer=get_earth(obstime))
hgs_frame = HeliographicStonyhurst(obstime=obstime)
hgc_frame = HeliographicCarrington(observer=get_earth(obstime))
hcc_frame = Heliocentric(observer=get_earth(obstime))


def get_point_chaincode(coord):
    if coord.frame.name == "heliocentric":
        components = coord.data.components[:-1]
    else:
        components = coord.data.components
    data = [copy.copy(getattr(coord.data, c).value) for c in components]
    if coord.frame.name == "heliocentric":
        data[-1] -= 90.0
    return f"POINT({' '.join([str(d) for d in data])})"


@pytest.mark.parametrize(("colname", "coord"), [
    ("hpc_coord", SkyCoord(*np.random.rand(2,len(obstime))*u.arcsec, frame=hpc_frame)),
    ("hgs_coord", SkyCoord(*np.random.rand(2,len(obstime))*u.deg, frame=hgs_frame)),
    ("hgc_coord", SkyCoord(*np.random.rand(2,len(obstime))*u.deg, frame=hgc_frame)),
    ("hrc_coord", SkyCoord(psi=np.random.rand(len(obstime))*360*u.deg,
                           rho=np.random.rand(len(obstime))*u.R_sun,
                           z=1*u.R_sun,
                           frame=hcc_frame,
                           representation_type='cylindrical')),
])
def test_chaincode_point_column_mapping(colname, coord):
    table = astropy.table.Table({
        "event_starttime": obstime.iso,
        colname: [get_point_chaincode(c) for c in coord],
    })
    _map_chain_code_columns_to_coordinates(table)
    assert isinstance(table[colname], SkyCoord)
    assert table[colname].frame.is_equivalent_frame(coord.frame)
    assert_quantity_allclose(table[colname].cartesian.xyz, coord.cartesian.xyz)


def get_polygon_chaincode(coord):
    if coord is None or coord == '':
        return coord
    c1 = copy.copy(getattr(coord.data, coord.data.components[0]).value)
    c2 = copy.copy(getattr(coord.data, coord.data.components[1]).value)
    if coord.frame.name == 'heliocentric':
        c2 -= 90.0
    data = [f'{c1} {c2}' for c1, c2 in zip(c1, c2)]
    return f"POLYGON(({','.join(data)}))"


@pytest.mark.parametrize(("colname", "coord"), [
    ("hpc_bbox", SkyCoord(*np.random.rand(2,7,len(obstime))*u.arcsec, frame=hpc_frame).T),
    ("hrc_bbox", SkyCoord(psi=np.random.rand(9,len(obstime))*360*u.deg,
                          rho=np.random.rand(9,len(obstime))*u.R_sun,
                          z=1*u.R_sun,
                          frame=hcc_frame,
                          representation_type='cylindrical').T),
    ("bound_chaincode", [SkyCoord(*np.random.rand(2,nb), unit='deg', frame='icrs') for nb in np.random.randint(6,11, size=len(obstime))]),
])
def test_chaincode_polygon_column_mapping(colname, coord):
    table = astropy.table.Table({
        "event_starttime": obstime.iso,
        "event_coordunit": len(obstime)*["deg"],
        colname: [get_polygon_chaincode(c) for c in coord],
    })
    _map_chain_code_columns_to_coordinates(table)
    if isinstance(table[colname], SkyCoord):
        assert table[colname].frame.is_equivalent_frame(coord.frame)
        assert_quantity_allclose(table[colname].cartesian.xyz, coord.cartesian.xyz)
    else:
        for i, c in enumerate(coord):
            assert isinstance(table[colname][i], SkyCoord)
            assert table[colname][i].frame.is_equivalent_frame(c.frame)
            assert_quantity_allclose(table[colname][i].cartesian.xyz, c.cartesian.xyz)
