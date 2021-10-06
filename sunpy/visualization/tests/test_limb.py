from matplotlib.figure import Figure

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.map.tests.conftest import aia171_test_map  # NoQA
from sunpy.visualization import draw_limb


def test_limb_invisible(aia171_test_map):
    aia_obs = aia171_test_map.observer_coordinate
    # Create a new observer on the opposite side of the Sun
    new_obs = SkyCoord(lon=aia_obs.lon + 180*u.deg,
                       lat=-aia_obs.lat,
                       radius=aia_obs.radius / 10,
                       frame=aia_obs.replicate_without_data())

    ax = Figure().add_subplot(111, projection=aia171_test_map)
    # Original observer
    visible, hidden = draw_limb(ax, aia_obs)
    assert visible is not None
    assert hidden is None
    # Far side observer
    visible, hidden = draw_limb(ax, new_obs)
    assert visible is None
    assert hidden is not None
