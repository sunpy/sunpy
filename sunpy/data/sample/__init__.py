"""SunPy sample data files"""
from __future__ import absolute_import

import os.path
from sunpy.util.net import url_exists
from os import rename, remove
from astropy.utils.data import download_file
from zipfile import ZipFile

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"

from sunpy import config
default_dir = config.get("downloads", "sample_dir")

_base_urls = ('http://data.sunpy.org/sample-data/', 'http://hesperia.gsfc.nasa.gov/~schriste/sunpy-sample-data/')

_files = {"AIA_171_IMAGE": "AIA20110319_105400_0171.fits",
          "RHESSI_IMAGE": "hsi_image_20101016_191218.fits",
          "EIT_195_IMAGE": "eit_l1_20020625_100011.fits",
          "CALLISTO_IMAGE": "BIR_20110922_103000_01.fit",
          "RHESSI_EVENT_LIST": "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
          "SWAP_LEVEL1_IMAGE": "swap_lv1_20120101_001607.fits",
          "AIA_193_IMAGE": "aia.lev1.193A_2013-09-21T16_00_06.84Z.image_lev1.fits.zip"
}

sample_files = {}
for key in _files:
    sample_files[key] = os.path.abspath(os.path.join(default_dir, _files[key]))

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
        for file_name in _files.itervalues():
            if url_exists(os.path.join(base_url, file_name)):
                f = download_file(os.path.join(base_url, file_name))
                
                real_name, ext = os.path.splitext(file_name)
                
                if ext == '.zip':
                    print("Unpacking: %s" % real_name)
                    with ZipFile(f, 'r') as zip_file:
                        zip_file.extract(real_name, default_dir)
                    remove(f)
                else:
                    # move files to the data directory
                    rename(f, os.path.join(default_dir, file_name))