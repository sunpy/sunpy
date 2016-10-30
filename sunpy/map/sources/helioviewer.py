"""Helioviewer source information.
"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1121

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

from astropy.visualization import LinearStretch

def source_stretch(meta, fits_stretch):
    """

    """
    # No helioviewer keyword in the map meta object
    if 'helioviewer' not in meta.keys():
        return fits_stretch
    elif meta['helioviewer']:
        # Helioviewer JPEG2000 files already have a stretched data values, so
        # just use a linear stretch.
        return LinearStretch()
    else:
        # Not a Helioviewer JPEG2000 file, so assume the data has not been
        # stretched and so use the FITS stretching as defined in the instrument
        # source.
        return fits_stretch
