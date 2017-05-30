# -*- coding: utf-8 -*-
"""SunPy sample data files"""
from __future__ import absolute_import, division, print_function

import os.path
import socket
import warnings
from zipfile import ZipFile
from shutil import move

from astropy.utils.data import download_file

from sunpy.extern import six
from sunpy.extern.six.moves.urllib.error import URLError

from sunpy.util.net import url_exists
from sunpy.util.config import get_and_create_sample_dir
from sunpy import config

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"


sampledata_dir = config.get("downloads", "sample_dir")

# urls to search for the sample data
_base_urls = (
    'http://data.sunpy.org/sample-data/',
    'https://github.com/sunpy/sunpy-sample-data/raw/master/')

# keys are file shortcuts
# values consist of filename as well as optional file extension if files are
# hosted compressed. This extension is removed after download.
_files = {
    "AIA_171_IMAGE": ("AIA20110319_105400_0171.fits", ""),
    "RHESSI_IMAGE": ("hsi_image_20101016_191218.fits", ""),
    "EIT_195_IMAGE": ("eit_l1_20020625_100011.fits", ""),
    "CALLISTO_IMAGE": ("BIR_20110922_103000_01.fit", ""),
    "RHESSI_EVENT_LIST": ("hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits", ""),
    "SWAP_LEVEL1_IMAGE": ("swap_lv1_20120101_001607.fits", ""),
    "AIA_193_IMAGE": ("aia.lev1.193A_2013-09-21T16_00_06.84Z.image_lev1.fits", ".zip"),
    "AIA_171_ROLL_IMAGE": ("aiacalibim5.fits.gz", ""),
    "AIA_94_CUTOUT": ("ssw_cutout_20121030_153001_AIA_94_.fts", ""),
    "EVE_LIGHTCURVE": ("20120620_EVE_L0CS_DIODES_1m.txt", ""),
    # Uncomment this if it needs to be used. Commented out to save bandwidth.
    # "LYRA_LIGHTCURVE": ("lyra_20110810-000000_lev2_std.fits.gz", ""),
    "LYRA_LEVEL3_LIGHTCURVE": ("lyra_20150101-000000_lev3_std.fits.gz", ""),
    "GOES_LIGHTCURVE": ("go1520120601.fits.gz", ""),
    "GBM_LIGHTCURVE": ("glg_cspec_n5_110607_v00.pha", ""),
    "NOAAINDICES_LIGHTCURVE": ("RecentIndices.txt", ""),
    "NOAAPREDICT_LIGHTCURVE": ("predicted-sunspot-radio-flux.txt", ""),
    "RHESSI_LIGHTCURVE": ("hsi_obssumm_20120601_018.fits.gz", ""),
    "NORH_LIGHTCURVE": ("tca110810", "")
}

sample_files = {}
for key in _files:
    sample_files[key] = os.path.abspath(os.path.join(sampledata_dir, _files[key][0]))

# Creating the directory for sample files to be downloaded
sampledata_dir = get_and_create_sample_dir()

def download_sample_data(progress=True, overwrite=False, timeout=None):
    """
    Download all sample data.
    
    Parameters
    ----------
    progress: `bool`
        Show a progress bar during download
    overwrite: `bool`
        If true, overwrites existing files.
    timeout: `float`
        The timeout in seconds. If `None` the default timeout is used from
        `astropy.utils.data.Conf.remote_timeout`.

    Returns
    -------
    None
    """
    print("Downloading all sample files to {}".format(sampledata_dir))
    for file_name in six.itervalues(_files):
        download_sample_file(file_name, url_list=_base_urls)


def download_sample_file(filename, url_list, progress=True, overwrite=False, timeout=None):
    """
    Download a sample data file and move it to the sample data directory. 
    Also, uncompresses the file if necessary.

    Parameters
    ----------
    filename: str
        Name of the file
    url_list: str or list
        urls where to look for the file
    progress: `bool`
        Show a progress bar during download
    overwrite: `bool`
        If exist overwrites the downloaded sample data.
    timeout: `float`
        The timeout in seconds. If `None` the default timeout is used from
        `astropy.utils.data.Conf.remote_timeout`.

    Returns
    -------
    None
    """

    if filename[-3:] == 'zip':
        uncompressed_filename = filename[:-4]
    else:
        uncompressed_filename = filename
    #print(uncompressed_filename)
    # check if the (uncompressed) file exists
    # print("Downloading sample files to {}".format(sampledata_dir))
    if not overwrite and os.path.isfile(os.path.join(sampledata_dir, uncompressed_filename)):
        print("File {} found and overwrite flag not set so skipping.".format(uncompressed_filename))
    else:
        # check each provided url to find the file
        for base_url in url_list:
            if base_url.count('github'):
                filename += '?raw=true'
            try:
                exists = url_exists(os.path.join(base_url, filename))
                if exists:
                    f = download_file(os.path.join(base_url, filename))
                    real_name, ext = os.path.splitext(f)

                    if ext == '.zip':
                        print("Unpacking: {}".format(real_name))
                        with ZipFile(f, 'r') as zip_file:
                            unzipped_f = zip_file.extract(real_name, sampledata_dir)
                        os.remove(f)
                        move(unzipped_f, os.path.join(sampledata_dir, uncompressed_filename))
                    else:
                        # move files to the data directory
                        move(f, os.path.join(sampledata_dir, uncompressed_filename))
                    # increment the number of files obtained to check later
                    break
            except (socket.error, socket.timeout) as e:
                warnings.warn("Download failed with error {}. \n"
                              "Retrying with different mirror.".format(e))
