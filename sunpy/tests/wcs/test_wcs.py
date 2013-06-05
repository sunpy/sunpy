from __future__ import absolute_import

#pylint: disable=E1103
import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_allclose

import sunpy
import sunpy.wcs as wcs

img = sunpy.make_map(sunpy.AIA_171_IMAGE)

# the following known_answers come from equivalent queries to IDL
# WCS implementation (http://hesperia.gsfc.nasa.gov/ssw/gen/idl/wcs/)

wcs.wcs.rsun_meters = img.rsun_meters

def test_convert_angle_units():
    actual = np.array([wcs._convert_angle_units(), wcs._convert_angle_units('arcsec'),
        wcs._convert_angle_units('arcmin'), wcs._convert_angle_units('deg'), 
        wcs._convert_angle_units('mas')])
    desired = np.array([np.deg2rad(1) / (60 * 60), np.deg2rad(1) / (60 * 60), 
        np.deg2rad(1) / 60.0, np.deg2rad(1), np.deg2rad(1) / (60 * 60 * 1000)])
    assert_allclose(actual, desired, rtol=1e-2, atol=0)
    
def test_conv_hpc_hcc():
    coord = [40.0, 32.0]
    result = wcs.convert_hpc_hcc(coord[0], coord[1], angle_units=img.units['x'])
    known_answer = [28748691, 22998953]
    assert_allclose(result, known_answer, rtol=1e-2, atol=0)
 
def test_conv_hcc_hpc():
    coord = [28748691, 22998953]
    result = wcs.convert_hcc_hpc(coord[0], coord[1], dsun_meters=img.dsun, 
        angle_units=img.units['x'])
    known_answer = [40.0, 32.0]
    assert_allclose(result, known_answer, rtol=1e-4, atol=0)
 
def test_conv_hcc_hg():
    coord = [13.0, 58.0]
    result = wcs.convert_hcc_hg(coord[0], coord[1], b0=img.heliographic_latitude, l0=img.heliographic_longitude)
    known_answer = [1.0791282e-06, -7.0640732]
    assert_allclose(result, known_answer, rtol=1e-1, atol=0)
 
def test_conv_hg_hcc():
    coord = [34.0, 96.0]
    result = wcs.convert_hg_hcc(coord[0], coord[1], b0=img.heliographic_latitude, 
                                l0=img.heliographic_longitude)
    known_answer = [-40653538.0, 6.7903529e8]
    assert_allclose(result, known_answer, rtol=1e-1, atol=0)
    
def test_conv_hg_hpc():
    coord = [34.0, 45.0]
    result = wcs.convert_hg_hpc(coord[0], coord[1], dsun_meters=img.dsun, 
                                b0=img.heliographic_latitude,
                                l0=img.heliographic_longitude, angle_units = img.units['x'])
    known_answer = [381.737592, 747.072612]
    assert_allclose(result, known_answer, rtol=1e-1, atol=0)
  
def test_conv_hpc_hg():
    coord = [381.737592, 747.072612]
    known_answer = [34.0, 45.0]
    result = wcs.convert_hpc_hg(coord[0], coord[1], dsun_meters=img.dsun, 
                                b0=img.heliographic_latitude,
                                l0=img.heliographic_longitude, angle_units = img.units['x'])
    assert_allclose(result, known_answer, rtol=1e-1, atol=0)
  