# -*- coding: utf-8 -*-
"""SunPy sample data files"""
from __future__ import absolute_import

from os import remove
import os.path
from zipfile import ZipFile
from urllib2 import URLError
from shutil import move

from astropy.utils.data import download_file

from sunpy.util.net import url_exists
from sunpy import config

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"


sampledata_dir = config.get("downloads", "sample_dir")

# urls to search for the sample data
_base_urls = (
    'http://data.sunpy.org/sample-data/',
    'http://hesperia.gsfc.nasa.gov/~schriste/sunpy-sample-data/',
    'https://github.com/ehsteve/sunpy-sample-data/raw/master/')

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
    "AIA_94_CUTOUT": ("ssw_cutout_20121030_153001_AIA_94_.fts", "")
}

sample_files = {}
for key in _files:
    sample_files[key] = os.path.abspath(os.path.join(sampledata_dir, _files[key][0]))


def download_sample_data(progress=True, overwrite=True):
    """
    Download the sample data.

    Parameters
    ----------
    progress: bool
        Show a progress bar during download
    overwrite: bool
        If exist overwrites the downloaded sample data.

    Returns
    -------
    None
    """
    number_of_files_fetched = 0
    print("Downloading sample files to " + sampledata_dir)
    for file_name in _files.itervalues():
        if not overwrite:
            if os.path.isfile(os.path.join(sampledata_dir,
                                           file_name[0])):
                number_of_files_fetched += 1
                continue

        for base_url in _base_urls:
            full_file_name = file_name[0] + file_name[1]
            if url_exists(os.path.join(base_url, full_file_name)):
                f = download_file(os.path.join(base_url, full_file_name))
                real_name, ext = os.path.splitext(full_file_name)

                if file_name[1] == '.zip':
                    print("Unpacking: %s" % real_name)
                    with ZipFile(f, 'r') as zip_file:
                        zip_file.extract(real_name, sampledata_dir)
                    remove(f)
                else:
                    # move files to the data directory
                    move(f, os.path.join(sampledata_dir, file_name[0]))
                # increment the number of files obtained to check later
                number_of_files_fetched += 1
                break

    if number_of_files_fetched < len(_files.keys()):
        raise URLError("Could not download all samples files. Problem with accessing sample data servers.")
