"""SunPy sample data files"""
from __future__ import absolute_import
import sunpy
import os
from zipfile import ZipFile
from os import remove
from os.path import exists, expanduser, isdir, join, splitext

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "sample")

def download(progress=True):
    """
    Download the sample data.
    """
    
    base_url = 'https://s3.amazonaws.com/bokeh_data/'
    
    files = [
        "AIA20110319_105400_0171.fits",
        "hsi_image_20101016_191218.fits",
        "eit_l1_20020625_100011.fits",
        "BIR_20110922_103000_01.fit",
        "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
        "swap_lv1_20120101_001607.fits"
        ]

    for file_name in files:
        _getfile(base_url, file_name, data_dir, progress=progress)

def _getfile(base_url, file_name, data_dir, progress=True):
    file_url = join(base_url, file_name)
    file_path = join(data_dir, file_name)

    url = urlopen(file_url)

    with open(file_path, 'wb') as file:
        file_size = int(url.headers["Content-Length"])
        print("Downloading: %s (%d bytes)" % (file_name, file_size))

        fetch_size = 0
        block_size = 16384

        while True:
            data = url.read(block_size)
            if not data:
                break

            fetch_size += len(data)
            file.write(data)

            if progress:
                status = "\r%10d [%6.2f%%]" % (fetch_size, fetch_size*100.0/file_size)
                stdout.write(status)
                stdout.flush()

    if progress:
        print()

    real_name, ext = splitext(file_name)

    if ext == '.zip':
        if not splitext(real_name)[1]:
            real_name += ".csv"

        print("Unpacking: %s" % real_name)

        with ZipFile(file_path, 'r') as zip_file:
            zip_file.extract(real_name, data_dir)

        remove(file_path)
#
# AIA20110319_105400_0171.fits
#
AIA_171_IMAGE = os.path.abspath(
    os.path.join(rootdir, "AIA20110319_105400_0171.fits")
)

#
# hsi_image_20101016_191218.fits
#
RHESSI_IMAGE = os.path.abspath(
    os.path.join(rootdir, "hsi_image_20101016_191218.fits")
)

#
# eit_l1_20020625_100011.fits
#
EIT_195_IMAGE = os.path.abspath(
    os.path.join(rootdir, "eit_l1_20020625_100011.fits")
)

CALLISTO_IMAGE = os.path.abspath(
    os.path.join(rootdir, "BIR_20110922_103000_01.fit")
)

# A stacked calibrated event list from RHESSI
# hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits
#
RHESSI_EVENT_LIST = os.path.abspath(
    os.path.join(rootdir, "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits")
)

#
# swap_lv1_20120101_001607.fits
#
SWAP_LEVEL1_IMAGE = os.path.join(rootdir, "swap_lv1_20120101_001607.fits")
