import pytest
import spiceypy

import astropy.units as u
from astropy.coordinates import ConvertError, SkyCoord, UnitSphericalRepresentation, frame_transform_graph
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import download_file

from sunpy.coordinates import spice
from sunpy.time import parse_time


@pytest.fixture(scope="module", params=[pytest.param(None, marks=pytest.mark.remote_data)])
def spice_test():
    # We use the cache for the SPK file because it's also used in other tests (does cache work?)
    spk = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
                        cache=True)

    # Leapseconds (LSK) and physical constants (PCK)
    lsk = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls")
    pck = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011_n0066.tpc")

    # Use frames (FK) from SunSPICE
    fk = download_file("https://soho.nascom.nasa.gov/solarsoft/packages/sunspice/data/heliospheric.tf")

    spice.initialize([spk, lsk, pck, fk])

    yield

    # The fact that we need to do this suggests that the cache is not persistent across tests...
    spiceypy.unload(spk)


def test_frame_creation(spice_test):
    expected = {"spice_ECLIPDATE", "spice_GSE", "spice_HCI", "spice_HEE", "spice_HEEQ", "spice_GEORTN"}
    assert expected.issubset(frame_transform_graph.get_names())


def test_double_initialize(spice_test):
    # Smoke test for an issue with calling initialize() more than once a session
    spice.initialize([])

    expected = {"spice_ECLIPDATE", "spice_GSE", "spice_HCI", "spice_HEE", "spice_HEEQ", "spice_GEORTN"}
    assert expected.issubset(frame_transform_graph.get_names())


def test_transformation(spice_test):
    coord = SkyCoord(1e7*u.km, 1e6*u.km, 1e5*u.km, representation_type='cartesian',
                     frame='spice_GSE', obstime='2023-10-17')
    new_coord = coord.spice_HEE
    new_coord.representation_type ='cartesian'

    assert_quantity_allclose(new_coord.x, 1.49128669e8*u.km - coord.x)
    assert_quantity_allclose(new_coord.y, -coord.y)
    assert_quantity_allclose(new_coord.z, coord.z)


def test_transformation_array_time(spice_test):
    coord = SkyCoord([1e7, 1e7]*u.km, [1e6, 1e6]*u.km, [1e5, 1e5]*u.km,
                     representation_type='cartesian',
                     frame='spice_GSE', obstime=['2023-10-17', '2023-10-18'])
    new_coord = coord.spice_HEE
    new_coord.representation_type ='cartesian'

    assert len(new_coord) == 2
    assert_quantity_allclose(new_coord.x, [1.49128669e8, 1.49085802e+08]*u.km - coord.x)
    assert_quantity_allclose(new_coord.y, -coord.y)
    assert_quantity_allclose(new_coord.z, coord.z)


def test_transformation_2d_okay(spice_test):
    coord = SkyCoord(30*u.deg, 45*u.deg, frame='spice_HCI', obstime='2023-10-17')
    new_coord = coord.spice_HEE

    assert isinstance(new_coord.data, UnitSphericalRepresentation)


def test_transformation_2d_not_okay(spice_test):
    coord = SkyCoord(30*u.deg, 45*u.deg, frame='spice_HCI', obstime='2023-10-17')
    with pytest.raises(ConvertError, match="due to a shift in origin"):
        assert coord.icrs

    coord = SkyCoord(30*u.deg, 45*u.deg, frame='icrs', obstime='2023-10-17')
    with pytest.raises(ConvertError, match="due to a shift in origin"):
        assert coord.spice_HCI


def test_get_body(spice_test):
    # Regression test
    earth_by_id = spice.get_body(399, '2023-10-17')
    earth_by_name = spice.get_body('earth', '2023-10-17')

    assert earth_by_id.name == 'icrs'
    assert_quantity_allclose(earth_by_id.ra, 21.3678774*u.deg)
    assert_quantity_allclose(earth_by_id.dec, 8.9883346*u.deg)
    assert_quantity_allclose(earth_by_id.distance, 1.47846987e8*u.km)

    assert earth_by_name == earth_by_id


def test_get_body_observer(spice_test):
    # Regression test
    earth = spice.get_body('earth', '2023-10-17')
    venus = spice.get_body('venus', '2023-10-17')
    venus_earth = spice.get_body('venus', '2023-10-17', observer=earth)
    venus_earth_hci = spice.get_body('venus', '2023-10-17', observer=earth, spice_frame_name='HCI')

    assert_quantity_allclose(venus_earth.separation_3d(venus), 11237.0641158*u.km)
    assert_quantity_allclose(venus_earth_hci.separation_3d(venus), 11237.0641158*u.km)


def test_get_body_array_time(spice_test):
    # Regression test
    obstime = parse_time(['2013-10-17', '2013-10-18'])
    earth = spice.get_body(399, obstime, spice_frame_name='HCI')

    assert len(earth) == 2
    assert earth[0].obstime == obstime[0]
    assert earth[1].obstime == obstime[1]
    assert earth[0].lon != earth[1].lon
    assert earth[0].lat != earth[1].lat
    assert earth[0].distance != earth[1].distance


def test_get_body_spice_frame(spice_test):
    # Regression test
    moon_gse = spice.get_body('moon', '2023-10-17', spice_frame_name='GSE')
    moon_gse.representation_type = 'cartesian'
    moon_hee = spice.get_body('moon', '2023-10-17', spice_frame_name='HEE')
    moon_hee.representation_type = 'cartesian'

    assert moon_gse.name == 'spice_GSE'
    assert_quantity_allclose(moon_gse.x, 349769.2488428*u.km)
    assert_quantity_allclose(moon_gse.y, 171312.08978192*u.km)
    assert_quantity_allclose(moon_gse.z, -14990.88338588*u.km)

    assert moon_hee.name == 'spice_HEE'
    assert_quantity_allclose(moon_hee.x, 1.487789e+08*u.km)
    assert_quantity_allclose(moon_hee.y, -171312.08978202*u.km)
    assert_quantity_allclose(moon_hee.z, -14990.88338588*u.km)
