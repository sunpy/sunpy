"""
This sub-module contains code for managing and downloading sample data
that's used only in sunpy documentation.
"""
from pathlib import Path

import pooch

from sunpy.map import Map
from sunpy.util.config import get_and_create_sample_dir

version = 'v1'
DOC_POOCH = pooch.create(
    # Use the default cache folder for the operating system
    path=Path(get_and_create_sample_dir()) / 'docs',
    # The remote data is on Github
    base_url="https://github.com/dstansby/sample-data/raw/venus-data/sunpy-docs/",
    version=version,
    registry={
        "aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1_lowres.fits": None,
    },
)


def fetch_venus_transit_map():
    """
    Load the AIA 1600 Venus transit map.
    """
    fname = DOC_POOCH.fetch("aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1_lowres.fits")
    return Map(fname)
