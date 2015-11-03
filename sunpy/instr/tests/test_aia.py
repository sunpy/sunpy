from __future__ import absolute_import

import tempfile

import numpy as np

import sunpy
import sunpy.data.test as test
from sunpy.instr.aia import aiaprep

# Define the original and prepped images first so they're available to all functions
original = sunpy.map.Map(test.aia_171_level1)
prep_map = aiaprep(original)


def test_aiaprep():
    # Test that header info for the map has been correctly updated
    # Check all of these for Map attributes and .meta values?
    # Check crpix values
    assert prep_map.meta['crpix1'] == prep_map.data.shape[1]/2.0 + 0.5
    assert prep_map.meta['crpix2'] == prep_map.data.shape[0]/2.0 + 0.5
    # Check cdelt values
    assert prep_map.meta['cdelt1']/0.6 == int(prep_map.meta['cdelt1']/0.6)
    assert prep_map.meta['cdelt2']/0.6 == int(prep_map.meta['cdelt2']/0.6)
    # Check rotation value, I am assuming that the inaccuracy in
    # the CROTA -> PCi_j matrix is causing the inaccuracy here
    np.testing.assert_allclose(prep_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert prep_map.meta['lvl_num'] == 1.5


def test_filesave():
    # Test that adjusted header values are still correct after saving the map
    # and reloading it.
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    prep_map.save(afilename, clobber=True)
    load_map = sunpy.map.Map(afilename)
    # Check crpix values
    assert load_map.meta['crpix1'] == prep_map.data.shape[1]/2.0 + 0.5
    assert load_map.meta['crpix2'] == prep_map.data.shape[0]/2.0 + 0.5
    # Check cdelt values
    assert load_map.meta['cdelt1']/0.6 == int(load_map.meta['cdelt1']/0.6)
    assert load_map.meta['cdelt2']/0.6 == int(load_map.meta['cdelt2']/0.6)
    # Check rotation value
    np.testing.assert_allclose(prep_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert load_map.meta['lvl_num'] == 1.5
