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

_base_urls = (
    'http://data.sunpy.org/sunpy/v1/',
    'https://github.com/sunpy/sample-data/raw/master/sunpy/v1/'
)

# Shortcut requirements:
# start with the instrument name then
# the wavelength or energy if needed then
# an optional description if needed then
# a reference name for the class into which the file will be opened
# (e.g. IMAGE for Maps, TIMESERIES for TimeSeries, SPECTRUM for Spectrum)
# All separated by underscores

# the files should include necessary extensions
_sample_files = {
    "AIA_131_IMAGE": "AIA20110607_063301_0131_lowres.fits",
    "AIA_171_IMAGE": "AIA20110607_063302_0171_lowres.fits",
    "AIA_211_IMAGE": "AIA20110607_063302_0211_lowres.fits",
    "AIA_335_IMAGE": "AIA20110607_063303_0335_lowres.fits",
    "AIA_094_IMAGE": "AIA20110607_063305_0094_lowres.fits",
    "AIA_1600_IMAGE": "AIA20110607_063305_1600_lowres.fits",
    "AIA_193_IMAGE": "AIA20110607_063307_0193_lowres.fits",
    "AIA_193_CUTOUT01_IMAGE": "AIA20110607_063307_0193_cutout.fits",
    "AIA_193_CUTOUT02_IMAGE": "AIA20110607_063931_0193_cutout.fits",
    "AIA_193_CUTOUT03_IMAGE": "AIA20110607_064555_0193_cutout.fits",
    "AIA_193_CUTOUT04_IMAGE": "AIA20110607_065219_0193_cutout.fits",
    "AIA_193_CUTOUT05_IMAGE": "AIA20110607_065843_0193_cutout.fits",
    "EIT_195_IMAGE": "eit_l1_20110607_203753.fits",
    "RHESSI_IMAGE": "hsi_image_20110607_063300.fits",
    "CALLISTO_SPECTRUM": "BIR_20110607_062400_10.fit",
    # Not in the sample-data repo
    # "RHESSI_EVENT_LIST": "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
    "SWAP_LEVEL1_IMAGE": "swap_lv1_20110607_063329.fits",
    "AIA_171_ROLL_IMAGE": "aiacalibim5.fits.gz",
    "EVE_TIMESERIES": "20110607_EVE_L0CS_DIODES_1m.txt",
    # Uncomment this if it needs to be used. Commented out to save bandwidth.
    # "LYRA_LIGHTCURVE": ("lyra_20110810-000000_lev2_std.fits.gz", ,
    "LYRA_LEVEL3_TIMESERIES": "lyra_20110607-000000_lev3_std.fits",
    "GOES_XRS_TIMESERIES": "go1520110607.fits",
    "GBM_TIMESERIES": "glg_cspec_n5_110607_v00.pha",
    "NOAAINDICES_TIMESERIES": "swpc_solar_cycle_indices.txt",
    "NOAAPREDICT_TIMESERIES": "predicted-sunspot-radio-flux.txt",
    "RHESSI_TIMESERIES": "hsi_obssumm_20110607_025.fits",
    "NORH_TIMESERIES": "tca110607.fits"
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
    for file_name in six.itervalues(_sample_files):
        get_sample_file(file_name, show_progress=show_progress,
                        url_list=_base_urls, overwrite=True)


def get_sample_file(filename, url_list, show_progress=True, overwrite=False,
                    timeout=None):
    """
    Downloads a sample file. Will download  a sample data file and move it to
    the sample data directory. Also, uncompresses zip files if necessary.
    Returns the local file if exists.

    Parameters
    ----------
    filename: `str`
        Name of the file
    url_list: `str` or `list`
        urls where to look for the file
    show_progress: `bool`
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
    if not overwrite and os.path.isfile(os.path.join(sampledata_dir,
                                                     uncompressed_filename)):
        return os.path.join(sampledata_dir, uncompressed_filename)
    else:
        # check each provided url to find the file
        for base_url in url_list:
            online_filename = filename
            if base_url.count('github'):
                online_filename += '?raw=true'
            try:
                url = six.moves.urllib_parse.urljoin(base_url, online_filename)
                exists = url_exists(url)
                if exists:
                    f = download_file(os.path.join(base_url, online_filename),
                                      show_progress=show_progress,
                                      timeout=timeout)
                    real_name, ext = os.path.splitext(f)

                    if ext == '.zip':
                        print("Unpacking: {}".format(real_name))
                        with ZipFile(f, 'r') as zip_file:
                            unzipped_f = zip_file.extract(real_name,
                                                          sampledata_dir)
                        os.remove(f)
                        move(unzipped_f, os.path.join(sampledata_dir,
                                                      uncompressed_filename))
                        return os.path.join(sampledata_dir,
                                            uncompressed_filename)
                    else:
                        # move files to the data directory
                        move(f, os.path.join(sampledata_dir,
                                             uncompressed_filename))
                        return os.path.join(sampledata_dir,
                                            uncompressed_filename)
            except (socket.error, socket.timeout) as e:
                warnings.warn("Download failed with error {}. \n"
                              "Retrying with different mirror.".format(e))
        # if reach here then file has not been downloaded.
        warnings.warn("File {} not found.".format(filename))
        return None
