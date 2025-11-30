import matplotlib.pyplot as plt
import pytest
from astropy.io import fits
from astropy.wcs import WCS

import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE

def test_crash_on_ape14_wcs_plotting():
    ndcube = pytest.importorskip('ndcube')
    from ndcube import NDCube

    image_data = fits.getdata(AIA_171_IMAGE)
    image_header = fits.getheader(AIA_171_IMAGE)

    cube = NDCube(image_data, WCS(image_header))[10:, 10:]
    aiamap = sunpy.map.Map(AIA_171_IMAGE)

    fig = plt.figure()
    ax = plt.subplot(projection=cube.wcs)
    
    aiamap.plot(axes=ax)
    plt.close(fig)