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

from sunpy.util.net import url_exists
from sunpy.util.config import get_and_create_sample_dir
from sunpy import config

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"


sampledata_dir = config.get("downloads", "sample_dir")

base_urls = (
    'http://data.sunpy.org/sample-data/',
    'https://github.com/sunpy/sunpy-sample-data/raw/master/'
)

# Shortcut requirements:
# the name of the class into which the file will opened must be included at the end of the file

# the files should include necessary extensions
files = {
    "AIA_171_IMAGE": "AIA20110319_105400_0171.fits",
    "RHESSI_IMAGE": "hsi_image_20101016_191218.fits",
    "EIT_195_IMAGE": "eit_l1_20020625_100011.fits",
    "CALLISTO_IMAGE": "BIR_20110922_103000_01.fit",
    "RHESSI_EVENT_LIST": "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
    "SWAP_LEVEL1_IMAGE": "swap_lv1_20120101_001607.fits",
    "AIA_193_IMAGE": "aia.lev1.193A_2013-09-21T16_00_06.84Z.image_lev1.fits.zip",
    "AIA_171_ROLL_IMAGE": "aiacalibim5.fits.gz",
    "AIA_94_CUTOUT_IMAGE": "ssw_cutout_20121030_153001_AIA_94_.fts",
    "EVE_TIMESERIES": "20120620_EVE_L0CS_DIODES_1m.txt",
    # Uncomment this if it needs to be used. Commented out to save bandwidth.
    # "LYRA_LIGHTCURVE": ("lyra_20110810-000000_lev2_std.fits.gz", ,
    "LYRA_LEVEL3_TIMESERIES": "lyra_20150101-000000_lev3_std.fits.gz",
    "GOES_XRS_TIMESERIES": "go1520120601.fits.gz",
    "GBM_TIMESERIES": "glg_cspec_n5_110607_v00.pha",
    "NOAAINDICES_TIMESERIES": "RecentIndices.txt",
    "NOAAPREDICT_TIMESERIES": "predicted-sunspot-radio-flux.txt",
    "RHESSI_TIMESERIES": "hsi_obssumm_20120601_018.fits.gz",
    "NORH_TIMESERIES": "tca110810"
}

# Creating the directory for sample files to be downloaded
sampledata_dir = get_and_create_sample_dir()


def download_sample_data(show_progress=True):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    show_progress: `bool`
        Show a progress bar during download

    Returns
    -------
    None
    """
    print("Downloading all sample files to {}. Overwriting if necessary.".format(sampledata_dir))
    for file_name in six.itervalues(files):
        get_sample_file(file_name, show_progress=show_progress, url_list=base_urls, overwrite=True)


def get_sample_file(filename, url_list, show_progress=True, overwrite=False, timeout=None):
    """
    Downloads a sample file. Will download  a sample data file and move it to the sample data directory.
    Also, uncompresses zip files if necessary. Returns the local file if exists.

    Parameters
    ----------
    filename: `str`
        Name of the file
    url_list: `str` or `list`
        urls where to look for the file
    progress: `bool`
        Show a progress bar during download
    overwrite: `bool`
        If True download and overwrite an existing file.
    timeout: `float`
        The timeout in seconds. If `None` the default timeout is used from
        `astropy.utils.data.Conf.remote_timeout`.

    Returns
    -------
    result: `str`
        The local path of the file. None if it failed.
    """

    if filename[-3:] == 'zip':
        uncompressed_filename = filename[:-4]
    else:
        uncompressed_filename = filename
    # check if the (uncompressed) file exists
    if not overwrite and os.path.isfile(os.path.join(sampledata_dir, uncompressed_filename)):
        return os.path.join(sampledata_dir, uncompressed_filename)
    else:
        # check each provided url to find the file
        for base_url in url_list:
            online_filename = filename
            if base_url.count('github'):
                online_filename += '?raw=true'
            try:
                exists = url_exists(os.path.join(base_url, online_filename))
                if exists:
                    f = download_file(os.path.join(base_url, online_filename), show_progress=show_progress)
                    real_name, ext = os.path.splitext(f)

                    if ext == '.zip':
                        print("Unpacking: {}".format(real_name))
                        with ZipFile(f, 'r') as zip_file:
                            unzipped_f = zip_file.extract(real_name, sampledata_dir)
                        os.remove(f)
                        move(unzipped_f, os.path.join(sampledata_dir, uncompressed_filename))
                        return os.path.join(sampledata_dir, uncompressed_filename)
                    else:
                        # move files to the data directory
                        move(f, os.path.join(sampledata_dir, uncompressed_filename))
                        return os.path.join(sampledata_dir, uncompressed_filename)
            except (socket.error, socket.timeout) as e:
                warnings.warn("Download failed with error {}. \n"
                              "Retrying with different mirror.".format(e))
        # if reach here then file has not been downloaded.
        warnings.warn("File {} not found.".format(filename))
        return None

