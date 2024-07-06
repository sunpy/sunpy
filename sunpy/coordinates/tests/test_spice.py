import numpy as np
import pytest
import spiceypy

import astropy.units as u
from astropy.coordinates import ConvertError, SkyCoord, UnitSphericalRepresentation, frame_transform_graph
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import download_file

from sunpy.coordinates import Helioprojective, get_earth, spice
from sunpy.time import parse_time


@pytest.fixture(scope="module", params=[pytest.param(None, marks=pytest.mark.remote_data)])
def spice_test():
    # We use the cache for the SPK file because it's also used in other tests (does cache work?)
    spk = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
                        cache=True)

    # Implicitly test adding a single kernel, but no kernel frames are installed here
    spice.initialize(spk)

    # Leapseconds (LSK) and physical constants (PCK)
    lsk = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls")
    pck = download_file("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011_n0066.tpc")

    # Use frames (FK) from SunSPICE
    fk = download_file("https://soho.nascom.nasa.gov/solarsoft/packages/sunspice/data/heliospheric.tf")

    # Use more frames (FK and IK) from Solar Orbiter
    solo_sc_fk = download_file("http://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/fk/solo_ANC_soc-sc-fk_V09.tf")
    solo_spd_ik = download_file("http://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/ik/solo_ANC_soc-epd-ik_V03.ti")

    # Implicitly test adding a kernel list
    spice.initialize([lsk, pck, fk, solo_sc_fk, solo_spd_ik])

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


def test_install_frame(spice_test):
    # Installing by ID is already tested through a successful initialize()

    assert "spice_IAU_SUN" not in frame_transform_graph.get_names()
    spice.install_frame('IAU_SUN')
    assert "spice_IAU_SUN" in frame_transform_graph.get_names()

    assert "spice_IAU_EARTH" not in frame_transform_graph.get_names()
    spice.install_frame('iau_earth')
    assert "spice_IAU_EARTH" in frame_transform_graph.get_names()


def test_install_frame_bad_id(spice_test):
    with pytest.raises(ValueError, match="not a valid SPICE frame ID"):
        spice.install_frame(999999)


def test_install_frame_bad_name(spice_test):
    with pytest.raises(ValueError, match="not a valid SPICE frame name"):
        spice.install_frame('NOT_A_FRAME_NAME')


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
    venus_earth_hci = spice.get_body('venus', '2023-10-17', observer=earth, spice_frame='HCI')

    assert_quantity_allclose(venus_earth.separation_3d(venus), 11237.0641158*u.km)
    assert_quantity_allclose(venus_earth_hci.separation_3d(venus), 11237.0641158*u.km)


def test_get_body_array_time(spice_test):
    # Regression test
    obstime = parse_time(['2013-10-17', '2013-10-18'])
    earth = spice.get_body(399, obstime, spice_frame='HCI')

    assert len(earth) == 2
    assert earth[0].obstime == obstime[0]
    assert earth[1].obstime == obstime[1]
    assert earth[0].lon != earth[1].lon
    assert earth[0].lat != earth[1].lat
    assert earth[0].distance != earth[1].distance


def test_get_body_spice_frame(spice_test):
    # Regression test
    moon_gse = spice.get_body('moon', '2023-10-17', spice_frame='GSE')
    moon_gse.representation_type = 'cartesian'
    moon_hee = spice.get_body('moon', '2023-10-17', spice_frame='HEE')
    moon_hee.representation_type = 'cartesian'

    assert moon_gse.name == 'spice_GSE'
    assert_quantity_allclose(moon_gse.x, 349769.2488428*u.km)
    assert_quantity_allclose(moon_gse.y, 171312.08978192*u.km)
    assert_quantity_allclose(moon_gse.z, -14990.88338588*u.km)

    assert moon_hee.name == 'spice_HEE'
    assert_quantity_allclose(moon_hee.x, 1.487789e+08*u.km)
    assert_quantity_allclose(moon_hee.y, -171312.08978202*u.km)
    assert_quantity_allclose(moon_hee.z, -14990.88338588*u.km)


def test_get_fov_rectangle(spice_test):
    fov = spice.get_fov('SOLO_EPD_STEP', '2023-10-17')
    width = 27*u.deg
    height = 13*u.deg
    # Because SPICE defines the rectangular FOV on a plane,
    # the height in spherical latitude is shorter at the corners
    adjusted_height = np.arctan(np.tan(height) * np.cos(width))
    assert_quantity_allclose(fov.lon, [1, -1, -1, 1] * width)
    assert_quantity_allclose(fov.lat, [1, 1, -1, -1] * adjusted_height)


def test_get_fov_circle(spice_test):
    obstime = parse_time(['2023-10-17', '2023-10-18'])
    fov = spice.get_fov('SOLO_EPD_SIS_ASW', obstime)
    angle = 11*u.deg
    assert fov.shape == (2, 100)
    assert (fov.obstime[:, 0] == obstime).all()
    assert_quantity_allclose(fov[0].cartesian.xyz, fov[1].cartesian.xyz)
    assert_quantity_allclose(fov[0, [0, 25, 50, 75]].lon, [1, 0, -1, 0] * angle, atol=1e-10*u.deg)
    assert_quantity_allclose(fov[0, [0, 25, 50, 75]].lat, [0, -1, 0, 1] * angle, atol=1e-10*u.deg)


def test_to_helioprojective(spice_test):
    # Define a coordinate in SunSPICE's GSE frame
    spice_coord = SkyCoord(UnitSphericalRepresentation(0.1*u.deg, 0.2*u.deg),
                           frame='spice_GSE', obstime='2023-10-17')
    hpc = spice_coord.to_helioprojective()

    assert isinstance(hpc.frame, Helioprojective)

    # The origin of spice_GSE (the Earth) should be used as the Helioprojective observer
    # We loosen tolerances because the planetary ephemeris is different
    earth = get_earth(hpc.obstime)
    assert_quantity_allclose(hpc.observer.lon, earth.lon, atol=1e-5*u.deg)
    assert_quantity_allclose(hpc.observer.lat, earth.lat, rtol=1e-6)
    assert_quantity_allclose(hpc.observer.radius, earth.radius, rtol=1e-6)

    # The angular offset from disk center should exactly match
    assert_quantity_allclose(np.arccos(hpc.cartesian.x), 0.22360671*u.deg)

    # We test the correctness of the HPC coordinate by transforming back to *our* GSE frame
    # and comparing against the original values in the SunSPICE GSE frame
    sunpy_coord = hpc.geocentricsolarecliptic
    assert_quantity_allclose(sunpy_coord.lon, spice_coord.lon, rtol=1e-6)
    assert_quantity_allclose(sunpy_coord.lat, spice_coord.lat, rtol=1e-6)

def test_get_rotation_matrix(spice_test):
    result1 = spice.get_rotation_matrix('HEEQ', 'J2000', '2024-07-04')
    expected_result1 = np.array([[ 0.20454304,  0.97118061,  0.12235349],
                                 [-0.87442558,  0.23746563, -0.42307208],
                                 [-0.43993415, -0.02045257,  0.8977971]])
    np.testing.assert_allclose(result1, expected_result1, atol=1e-6)

    result2 = spice.get_rotation_matrix('HEEQ', 'HEEQ', '2024-07-04', '2024-07-18')
    expected_result2 = np.array([[ 0.97314148,  0.23020788,  0],
                                 [-0.23020788,  0.97314148, 0],
                                 [0, 0, 1]])
    np.testing.assert_allclose(result2, expected_result2, atol=1e-6)
