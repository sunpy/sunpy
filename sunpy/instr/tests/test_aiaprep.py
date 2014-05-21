from __future__ import absolute_import

import sunpy
from sunpy.instr.aia import aiaprep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Define test image first so it's accessable to all functions.
#original = sunpy.map.Map(sunpy.AIA_171_IMAGE)
data = np.zeros((101, 101))
data[50:60, 50:80] = 1.0
# An entirely made up test header.
header = {'cdelt1': 0.6,
          'cdelt2': 0.6,
          #crota1?
          #crota2: 
          'crpix1': data.shape[1]/2.0 + 0.5,
          'crpix2': data.shape[0]/2.0 + 0.5,
          #crval1:
          #crval2:
          'cunit1': 'arcsec',
          'cunit2': 'arcsec',
          #inst_rot: 
          'lvl_num': 1.0,
          'naxis': 2,
          'naxis1': data.shape[1],
          'naxis2': data.shape[0], # Hope these are the right way round...
          'date-obs': '2014-05-20 14:35'}

original = sunpy.map.Map(data, header)
"""
fig = plt.figure()
original.plot(cmap=cm.binary)
plt.show()

header = original.meta.keys()
header.sort()
for item in header:
    print item, original.meta[item]
"""

# Tolerance for tests - high for now just so they run. I'll make it much smaller
# when actually testing
atol = 0.1

def plot_results(expect, result, diff):
    fig = plt.figure()
    
    fig.add_subplot(1, 3, 1)
    expect.plot(cmap=cm.binary)
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 2)
    result.plot(cmap=cm.binary)
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 3)
    diff.plot(cmap=cm.coolwarm)
    plt.colorbar(orientation='horizontal')


def test_aiaprep_rotation():
    # Check 360 degree rotation against original image
    original.meta['crota2'] = 360.0
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(original.data-rot_map.data, rot_map.meta)
    plot_results(original, rot_map, diff_map)
    assert np.allclose(original.data, rot_map.data, atol=atol)
    assert abs(original.data.mean() - rot_map.data.mean()) <= atol
    plt.close()
    
    # Check incremental 360 degree rotation against original image

    # Test 90 degree rotation against expected outcome
    original.meta['crota2'] = 90.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy()), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-rot_map.data, rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    assert np.allclose(expected.data, rot_map.data, atol=atol)
    assert abs(expected.data.mean() - rot_map.data.mean()) <= atol
    plt.close()

    # Test 90 degree rotation against -270 degree rotation
    original.meta['crota2'] = -270.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy()), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-rot_map.data, rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    assert np.allclose(expected.data, rot_map.data, atol=atol)
    assert abs(expected.data.mean() - rot_map.data.mean()) <= atol
    plt.close()

    # Test -90 degree rotation against 270 degree rotation
    original.meta['crota2'] = 90.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy(), -3), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-rot_map.data, rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    assert np.allclose(expected.data, rot_map.data, atol=atol)
    assert abs(expected.data.mean() - rot_map.data.mean()) <= atol
    plt.close()

    # Check rotated and derotated image against original
    original.meta['crota2'] = 10.0
    rot_map = aiaprep(original)
    rot_map.meta['crota2'] = -10.0
    rot_map = aiaprep(rot_map)
    diff_map = sunpy.map.GenericMap(original.data-rot_map.data, rot_map.meta)
    plot_results(original, rot_map, diff_map)
    assert np.allclose(original.data, rot_map.data, atol=atol)
    assert abs(original.data.mean() - rot_map.data.mean()) <= atol
    plt.close()


def test_aiaprep_shift():
    # Check a shifted shape against expected outcome
    original.meta['crota2'] = 0.0
    expect_data = np.zeros((101, 101))
    expect_data[45:55, 35:65] = 1.0 # Not entirely convinced by these numbers
    expected = sunpy.map.GenericMap(expect_data, original.meta.copy())
    original.meta['crpix1'], original.meta['crpix2'] = 65, 55
    shift_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-shift_map.data, shift_map.meta)
    plot_results(expected, shift_map, diff_map)
    assert np.allclose(expected.data, shift_map.data, atol=atol)
    assert abs(expected.data.mean() - shift_map.data.mean()) <= atol
    plt.close()

    # Check shifted and unshifted shape against original image
    shift_map = aiaprep(original)
    shift_map.meta['crpix1'], shift_map.meta['crpix2'] = data.shape[1]/2.0 + 0.5, data.shape[0]/2.0 + 0.5
    shift_map = aiaprep(shift_map)
    diff_map = sunpy.map.GenericMap(original.data-shift_map.data, shift_map.meta)
    plot_results(original, shift_map, diff_map)
    assert np.allclose(original.data, shift_map.data, atol=atol)
    assert abs(original.data.mean() - shift_map.data.mean()) <= atol
    plt.close()


def test_aiaprep_scale():
    # Check a scaled image against the expected outcome
    original.meta['crota2'] = 0.0
    expect_data = np.zeros((101, 101))
    #expect_data[??, ??] = 1.0 # Not entirely sure what the expected outcome should be
    expected = sunpy.map.GenericMap(expect_data, original.meta.copy())
    original.meta['cdelt1'], original.meta['cdelt2'] = 0.3, 0.3
    scale_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-scale_map.data, scale_map.meta)
    plot_results(expected, scale_map, diff_map)
    assert np.allclose(expected.data, scale_map.data, atol=atol)
    assert abs(expected.data.mean() - scale_map.data.mean()) <= atol
    plt.close()
    
    # Check a scaled and descaled image against the original
    scale_map = aiaprep(original)
    scale_map.meta['cdelt1'], scale_map.meta['cdelt2'] = 1.2, 1.2 # I think this should be the same as changing the target scale in aiaprep()
    scale_map = aiaprep(scale_map)
    diff_map = sunpy.map.GenericMap(original.data-scale_map.data, scale_map.meta)
    plot_results(original, scale_map, diff_map)
    assert np.allclose(original.data, scale_map.data, atol=atol)
    assert abs(original.data.mean() - scale_map.data.mean()) <= atol
    plt.close()


def test_aiaprep_all():
    # Check a shifted, rotated and scaled shape against expected outcome
    original.meta['crota2'] = 90.0
    expect_data = np.zeros((101, 101))
    #expect_data[??, ??] = 1.0 # Not entirely sure what the expected outcome should be
    expected = sunpy.map.GenericMap(expect_data, original.meta.copy())
    original.meta['crpix1'], original.meta['crpix2'] = 65, 55 # Not entirely sure about this either
    original.meta['cdelt1'], original.meta['cdelt2'] = 0.3, 0.3
    prep_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(expected.data-prep_map.data, prep_map.meta)
    plot_results(expected, prep_map, diff_map)
    assert np.allclose(expected.data, prep_map.data, atol=atol)
    assert abs(expected.data.mean() - prep_map.data.mean()) <= atol
    plt.close()

    # Check a prepped and de-prepped shape against original image
    original.meta['crota2'] = 10.0
    prep_map = aiaprep(original)
    prep_map.meta['crota2'] = -10.0
    prep_map.meta['crpix1'], original.meta['crpix2'] = data.shape[1]/2.0 + 0.5, data.shape[0]/2.0 + 0.5 # Not entirely sure about this either
    prep_map.meta['cdelt1'], original.meta['cdelt2'] = 1.2, 1.2
    prep_map = aiaprep(prep_map)
    diff_map = sunpy.map.GenericMap(original.data-prep_map.data, prep_map.meta)
    plot_results(original, prep_map, diff_map)
    assert np.allclose(expected.data, prep_map.data, atol=atol)
    assert abs(expected.data.mean() - prep_map.data.mean()) <= atol
    plt.close()

    # Test that header info for the map has been correctly updated
    # Check all of these for Map attributes and .meta values?
    # Also may be worth checking they stay the same when saved, I'm sure I've had issues with that before.
    # Check crpix values
    assert prep_map.meta['crpix1'] == prep_map.shape[1]/2.0 + 0.5
    assert prep_map.meta['crpix2'] == prep_map.shape[1]/2.0 + 0.5
    # Check cdelt values
    assert prep_map.meta['cdelt1'] == 0.6
    assert prep_map.meta['cdelt2'] == 0.6
    # Check crota values
    assert prep_map.meta['crota1'] == 0.0
    assert prep_map.meta['crota2'] == 0.0
    # Check level number
    assert prep_map.meta['lvl_num'] == 1.5


try:
    test_aiaprep_rotation()
    test_aiaprep_shift()
    test_aiaprep_scale()
    test_aiaprep_all()
except AssertionError:
    plt.show()
    raise