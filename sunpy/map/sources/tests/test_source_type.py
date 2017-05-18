"""
Test cases distinguishing the source types.
"""

import os
import glob
from astropy.visualization import LinearStretch
from sunpy.map.sources.source_type import from_helioviewer_project, source_stretch
from sunpy.map import Map
import sunpy.data.test

path = sunpy.data.test.rootdir
fitspath = glob.glob(os.path.join(path, "aia_171_level1.fits"))
aia = Map(fitspath)

fitspath = glob.glob(os.path.join(path, "2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2"))
hvjp2 = Map(fitspath)


def test_from_helioviewer_project():
    """Tests """
    assert not from_helioviewer_project(aia.meta)
    assert from_helioviewer_project(hvjp2.meta)


def test_source_stretch():
    aia_fits_stretch = aia.plot_settings['norm'].stretch
    assert source_stretch(aia.meta, aia_fits_stretch) is aia_fits_stretch
    assert isinstance(source_stretch(hvjp2.meta, aia_fits_stretch), LinearStretch)

