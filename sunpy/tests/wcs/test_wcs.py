from __future__ import absolute_import

#pylint: disable=E1103
import numpy as np
from numpy.testing import assert_array_almost_equal
import pyfits

import sunpy
import sunpy.wcs as wcs

fits = pyfits.open(sunpy.AIA_171_IMAGE)
header = fits[0].header
img = sunpy.make_map(sunpy.AIA_171_IMAGE)

# the following known_answers come from equivalent queries to IDL
# WCS implementation (http://hesperia.gsfc.nasa.gov/ssw/gen/idl/wcs/)

def test_conv_hpc_hcc():
    coord = [40.0, 32.0]
    result = wcs.convert_hpc_hcc(img.rsun_meters, 
                                 img.dsun, img.units['x'], img.units['y'], 
                                 coord[0], coord[1])
    known_answer = [28748691, 22998953]
    magnitude = np.floor(np.log10(np.abs(known_answer)))
    assert_array_almost_equal(result*10**(-magnitude), 
                              known_answer*10**(-magnitude), decimal=2)
 
def test_conv_hcc_hpc():
    coord = [34.0, 132.0]
    result = wcs.convert_hcc_hpc(img.rsun_meters, img.dsun, 
                                 coord[0], coord[1])
    known_answer = [1.3140782e-08, 5.1017152e-08]
    magnitude = np.floor(np.log10(np.abs(known_answer)))
    assert_array_almost_equal(result*10**(-magnitude), 
                              known_answer*10**(-magnitude), decimal=2)
 
def test_conv_hcc_hg():
    coord = [13.0, 58.0]
    result = wcs.convert_hcc_hg(img.rsun_meters, 
                                img.heliographic_latitude, img.heliographic_longitude,
                                coord[0], coord[1])
    known_answer = [1.0791282e-06, -7.0640732]
    magnitude = np.floor(np.log10(np.abs(known_answer)))
    assert_array_almost_equal(result*10**(-magnitude), 
                              known_answer*10**(-magnitude), decimal=2)
 

def test_conv_hg_hcc():
    coord = [34.0, 96.0]
    result = wcs.convert_hg_hcc(img.rsun_meters, img.heliographic_latitude, 
                                img.heliographic_longitude, coord[0], coord[1])
    known_answer = [-40653538.0, 6.7903529e08]
    magnitude = np.floor(np.log10(np.abs(known_answer)))
    assert_array_almost_equal(result*10**(-magnitude), 
                              known_answer*10**(-magnitude), decimal=2)
    
def test_conv_hg_hpc():
    coord = [34.0, 45.0]
    result = wcs.convert_hg_hpc(img.rsun_meters, img.dsun, 
                                img.heliographic_latitude,
                                img.heliographic_longitude,
                                coord[0], coord[1])
    known_answer = [0.10603822, 0.20752017]
    magnitude = np.floor(np.log10(known_answer))
    assert_array_almost_equal(result*10**(-magnitude), 
                              known_answer*10**(-magnitude), decimal=2)
  
