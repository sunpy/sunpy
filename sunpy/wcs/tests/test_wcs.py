from __future__ import absolute_import

#pylint: disable=E1103
import numpy as np
from numpy.testing import assert_allclose

import astropy.units as u

import sunpy
import sunpy.map
import sunpy.wcs as wcs
import sunpy.sun as sun

from numpy.testing import assert_array_almost_equal

img = sunpy.map.Map(sunpy.AIA_171_IMAGE)

# the following known_answers come from equivalent queries to IDL
# WCS implementation (http://hesperia.gsfc.nasa.gov/ssw/gen/idl/wcs/)

wcs.wcs.rsun_meters = img.rsun_meters

#def test_convert_pixel_to_data():
#    actual =
#    reference_pixel = [img.reference_pixel['x'], img.reference_pixel['y']]
#    reference_coordinate = [img.reference_coordinate['x'], img.reference_coordinate['y']]
#    scale = [img.scale['x'], img.scale['y']]
#    actual = wcs.convert_pixel_to_data(scale,img.shape, reference_pixel, reference_coordinate, 0, 0)


def test_conv_hpc_hcc():
    coord = [40.0 * u.arcsec, 32.0 * u.arcsec]
    result = wcs.convert_hpc_hcc(coord[0], coord[1])
    known_answer = [28748691 * u.meter, 22998953 * u.meter]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # Test the dsun_meters parameter for a distance of 0.5 AU
    dist = 0.5 * sun.constants.au.si
    result = wcs.convert_hpc_hcc(coord[0], coord[1], dist)
    known_answer = [14370494 * u.meter, 11496395 * u.meter]
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)

    # Make sure that z coordinate is returned if parameter z is True
    result = wcs.convert_hpc_hcc(coord[0], coord[1], z=True)
    known_answer = [28748691*u.meter, 22998953*u.meter, 695016924*u.meter]
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)

def test_conv_hcc_hpc():
    coord = [28748691*u.meter, 22998953*u.meter]
    result = wcs.convert_hcc_hpc(coord[0], coord[1], dsun_meters=img.dsun)
    known_answer = [40.0*u.arcsec.to(u.deg), 32.0*u.arcsec.to(u.deg)]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol = 0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol = 0)

def test_conv_hcc_hg():
    coord = [13.0*u.meter, 58.0*u.meter]
    result = wcs.convert_hcc_hg(coord[0], coord[1], b0_deg=img.heliographic_latitude, l0_deg=img.heliographic_longitude)
    known_answer = [1.0791282e-06 * u.deg, -7.0640732 * u.deg]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # Make sure that r value is returned if radius=True
    result = wcs.convert_hcc_hg(coord[0], coord[1], b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude, radius=True)
    known_answer = [1.0791282e-06 * u.deg, -7.0640732 * u.deg, sun.constants.radius.si]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

def test_conv_hg_hcc():
    coord = [34.0 * u.deg, 96.0 * u.deg]
    result = wcs.convert_hg_hcc(coord[0], coord[1], b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude)
    known_answer = [-40653538.0 * u.meter, 6.7903529e8 * u.meter]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # Test the radius parameter using half of the Sun's radius
    known_answer = [x / 2.0 for x in known_answer]
    radius = sun.constants.radius.si / 2.0
    result = wcs.convert_hg_hcc(coord[0], coord[1], b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude, r=radius)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # Make sure that z coordinates are returned if z=True
    known_answer = [-40653538.0 * u.meter, 6.7903529e8 * u.meter, -1.4487837e8 * u.meter]
    result = wcs.convert_hg_hcc(coord[0], coord[1], b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude, z=True)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # If z < 0, using occultation should make the return coordinates nan
    coord2 = [55.0 * u.deg, 56.0 * u.deg]
    known_answer = [[np.nan, 3.1858718e8 * u.meter], [np.nan, 5.9965928e8 * u.meter]]
    coords = zip(coord, coord2)
    result = wcs.convert_hg_hcc(*coords, b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude, occultation=True)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

def test_conv_hg_hpc():
    coord = [34.0 * u.deg, 45.0 * u.deg]
    result = wcs.convert_hg_hpc(coord[0], coord[1], dsun_meters=img.dsun,
                                b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude)
    known_answer = [381.737592 * u.arcsec.to(u.deg), 747.072612 * u.arcsec.to(u.deg)]
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

    # Test to make sure occultation parameter works
    coord = [34.0, 96.0]
    coord2 = [55.0, 56.0]
    known_answer = [[np.nan, 441.65710359], [np.nan, 831.30194808]]
    coords = zip(coord, coord2)
    result = wcs.convert_hg_hpc(*coords, dsun_meters=img.dsun,
                b0_deg=img.heliographic_latitude, l0_deg=img.heliographic_longitude,
                occultation=True)
    assert_allclose(result, known_answer, rtol=1e-2, atol=0)

def test_conv_hpc_hg():
    coord = [382 * u.arcsec, 748 * u.arcsec]
    known_answer = [34.091299 * u.deg, 45.095130 * u.deg]
    result = wcs.convert_hpc_hg(coord[0], coord[1], dsun_meters=img.dsun,
                                b0_deg=img.heliographic_latitude,
                                l0_deg=img.heliographic_longitude)
    assert_allclose(result[0], known_answer[0], rtol=1e-2, atol=0)
    assert_allclose(result[1], known_answer[1], rtol=1e-2, atol=0)

def test_convert_to_coord():
    x, y = (34.0, 96.0)
    b0_deg = img.heliographic_latitude
    l0_deg = img.heliographic_longitude
    dsun=img.dsun
    def check_conversion(from_coord, to_coord, expected):
        # Make sure that wcs.convert_to_coord returns the expected value
        assert_allclose(wcs.convert_to_coord(x * u.deg, y * u.deg, from_coord, to_coord,
            b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun)[0],
            expected[0], rtol=1e-2, atol=0)
        assert_allclose(wcs.convert_to_coord(x * u.deg, y * u.deg, from_coord, to_coord,
            b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun)[1],
            expected[1], rtol=1e-2, atol=0)        
    check_conversion('hpc', 'hg', wcs.convert_hpc_hg(x*u.deg, y*u.deg, b0_deg=b0_deg,
                        l0_deg=l0_deg, dsun_meters=dsun))
    check_conversion('hg', 'hcc', wcs.convert_hg_hcc(x*u.deg, y*u.deg,
                        b0_deg=b0_deg, l0_deg=l0_deg))
    check_conversion('hg', 'hpc', wcs.convert_hg_hpc(x*u.deg, y*u.deg, b0_deg=b0_deg,
                        l0_deg=l0_deg, dsun_meters=dsun))
    check_conversion('hpc', 'hcc', wcs.convert_hpc_hcc(x*u.deg, y*u.deg, dsun_meters=dsun))

    def check_conversion(from_coord, to_coord, expected):
        assert_allclose(wcs.convert_to_coord(x * u.meter, y * u.meter, from_coord, to_coord,
            b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun)[0],
            expected[0], rtol=1e-2, atol=0)
        assert_allclose(wcs.convert_to_coord(x * u.meter, y * u.meter, from_coord, to_coord,
            b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun)[1],
            expected[1], rtol=1e-2, atol=0) 

    check_conversion('hcc', 'hg', wcs.convert_hcc_hg(x*u.meter, y*u.meter,
                        b0_deg=b0_deg, l0_deg=l0_deg))
    check_conversion('hcc', 'hpc', wcs.convert_hcc_hpc(x * u.meter,
                        y * u.meter, dsun_meters=dsun))   
    

def test_convert_back():
    # Make sure transformation followed by inverse transformation returns
    # the original coordinates
    coord = [40.0 * u.deg, 32.0 * u.deg]
    assert_allclose(wcs.convert_hcc_hpc(*wcs.convert_hpc_hcc(*coord))[0],
                    coord[0], rtol=1e-2, atol=0)
    assert_allclose(wcs.convert_hcc_hpc(*wcs.convert_hpc_hcc(*coord))[1],
                    coord[1], rtol=1e-2, atol=0)    
    coord = [13.0 * u.meter, 58.0 * u.meter]
    assert_allclose(wcs.convert_hg_hcc(*wcs.convert_hcc_hg(*coord)),
                    coord, rtol=1e-2, atol=0)
    coord = [34.0 * u.deg, 45.0 * u.deg]
    assert_allclose(wcs.convert_hpc_hg(*wcs.convert_hg_hpc(*coord)),
                    coord, rtol=1e-2, atol=0)
