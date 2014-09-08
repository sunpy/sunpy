"""SunPy sample data files"""
from __future__ import absolute_import

import os
from os.path import join, splitext
from os import rename, remove
from astropy.utils.data import download_file
from zipfile import ZipFile
import urllib2

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"

from sunpy import config
default_dir = config.get("downloads", "sample_dir")

_base_urls = ('http://data.sunpy.org/sample-data/', 'http://hesperia.gsfc.nasa.gov/~schriste/sunpy-sample-data/')

_files = [
        "AIA20110319_105400_0171.fits",
        "hsi_image_20101016_191218.fits",
        "eit_l1_20020625_100011.fits",
        "BIR_20110922_103000_01.fit",
        "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
        "swap_lv1_20120101_001607.fits",
        "aia.lev1.193A_2013-09-21T16_00_06.84Z.image_lev1.fits.zip"
        ]

def download(progress=True):
    """
    Download the sample data.
    
    Parameters
    ----------
    progress: bool
        Show a progress bar during download
    
    Returns
    -------
    None
    """
    
    for base_url in _base_urls:
        for file_name in _files:
            if url_exists(join(base_url, file_name)):
                f = download_file(join(base_url, file_name))
                
                real_name, ext = splitext(file_name)
                
                if ext == '.zip':
                    print("Unpacking: %s" % real_name)
                    with ZipFile(f, 'r') as zip_file:
                        zip_file.extract(real_name, default_dir)
                    remove(f)
                else:
                    # move files to the data directory
                    rename(f, join(default_dir, file_name))
                    
def url_exists(url, timeout=2):
    """
    Checks whether a url is online.

    Parameters
    ----------
    url: str
        A string containing a URL

    Returns
    -------
    value: bool

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> url_exists('http://www.google.com')
    True
    >>> url_exists('http://aslkfjasdlfkjwerf.com')
    False
    """
    try:
        urllib2.urlopen(url, timeout=timeout)
    except urllib2.HTTPError, e:
        #print(url)
        #print(e.reason)
        return False
    except urllib2.URLError, e:
        #print(url)
        #print(e.reason)
        return False
    else:
        return True

#
# AIA20110319_105400_0171.fits
#
AIA_171_IMAGE = os.path.abspath(
    join(default_dir, "AIA20110319_105400_0171.fits")
)

#
# hsi_image_20101016_191218.fits
#
RHESSI_IMAGE = os.path.abspath(
    join(default_dir, "hsi_image_20101016_191218.fits")
)

#
# eit_l1_20020625_100011.fits
#
EIT_195_IMAGE = os.path.abspath(
    join(default_dir, "eit_l1_20020625_100011.fits")
)

CALLISTO_IMAGE = os.path.abspath(
    join(default_dir, "BIR_20110922_103000_01.fit")
)

# A stacked calibrated event list from RHESSI
# hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits
#
RHESSI_EVENT_LIST = os.path.abspath(
    join(default_dir, "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits")
)

#
# swap_lv1_20120101_001607.fits
#
SWAP_LEVEL1_IMAGE = join(default_dir, "swap_lv1_20120101_001607.fits")
