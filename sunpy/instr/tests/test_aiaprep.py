from __future__ import absolute_import

import sunpy
from sunpy.instr.aia import aiaprep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Define test image first so it's accessable to all functions.
#original = sunpy.map.Map(sunpy.AIA_171_IMAGE)
data = np.zeros((100, 100))
data[50:60, 50:80] = 1.0
# An entirely made up test header.
header = {'cdelt1': 0.599,
          'cdelt2': 0.599,
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

fig = plt.figure()
original.plot(cmap=cm.binary)
plt.show()

header = original.meta.keys()
header.sort()
for item in header:
    print item, original.meta[item]


def test_aiaprep_rotation():
    # Check 360 degree rotation against original image
    original.meta['crota2'] = 360.0
    rot_map = aiaprep(original)
    fig = plt.figure()
    fig.add_subplot(1, 2, 1)
    original.plot(cmap=cm.binary)
    fig.add_subplot(1, 2, 2)
    rot_map.plot(cmap=cm.binary)
    assert np.allclose(original.data, rot_map.data)
    # Check incremental 360 degree rotation against original image
    # Test 90 degree rotation against expected outcome
    # Test 90 degree rotation against -270 degree rotation
    # Test -90 degree rotation against 270 degree rotation
    # Check rotated and derotated image against original


def test_aiaprep_shift():
    # Check a shifted shape against expected outcome
    # Check shifted and unshifted shape against original image
    pass


def test_aiaprep_rot_and_shift():
    # Check a shifted and rotated shape against expected outcome
    # Check a shifted, rotated, then unshifted and derotated shape against original image
    pass


def test_aiaprep_headerinfo():
    # Test that header info for the map has been correctly updated
    pass

"""In addition to all of the above, I need to test the scaling as well.
Remembering that it hasn't neccessarily worked if the header says it has.
"""

try:
    test_aiaprep_rotation()
except AssertionError:
    plt.show()
    raise