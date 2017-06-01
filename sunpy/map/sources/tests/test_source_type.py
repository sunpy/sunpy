"""
Test cases distinguishing the source types.
"""

import os
import glob
from astropy.visualization import LinearStretch
from sunpy.map.sources.source_type import from_helioviewer_project, source_stretch
from sunpy.map import Map
import sunpy.data.test
from sunpy.tests.helpers import skip_glymur

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "aia_171_level1.fits"))
aia = Map(fitspath)
jp2path = glob.glob(os.path.join(path, "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))


@skip_glymur
def test_from_helioviewer_project():
    """Tests if we are able to determine if a file is from the Helioviewer
    Project or not."""
    hvjp2 = Map(jp2path)
    assert not from_helioviewer_project(aia.meta)
    assert from_helioviewer_project(hvjp2.meta)


@skip_glymur
def test_source_stretch():
    """
    Tests that the correct stretch function is returned.
    """
    hvjp2 = Map(jp2path)
    aia_fits_stretch = aia.plot_settings['norm'].stretch
    assert source_stretch(aia.meta, aia_fits_stretch) is aia_fits_stretch
    assert isinstance(source_stretch(hvjp2.meta, aia_fits_stretch), LinearStretch)

