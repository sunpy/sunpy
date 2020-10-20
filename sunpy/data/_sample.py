import os
from pathlib import Path
from urllib.parse import urljoin

from sunpy import log
from sunpy.util.config import _is_writable_dir, get_and_create_sample_dir
from sunpy.util.parfive_helpers import Downloader

_BASE_URLS = (
    'https://github.com/sunpy/sample-data/raw/master/sunpy/v1/',
    'http://data.sunpy.org/sunpy/v1/',
)

# Shortcut requirements:
# start with the instrument name then
# the wavelength or energy if needed then
# an optional description if needed then
# a reference name for the class into which the file will be opened
# (e.g. IMAGE for Maps, TIMESERIES for TimeSeries, SPECTRUM for Spectrum)
# All separated by underscores

# the files should include necessary extensions
_SAMPLE_DATA = {
    # Do roll image first because it's the largest file.
    "AIA_171_ROLL_IMAGE": "aiacalibim5.fits.gz",
    "HMI_LOS_IMAGE": "HMI20110607_063211_los_lowres.fits",
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
    "SWAP_LEVEL1_IMAGE": "swap_lv1_20110607_063329.fits",
    "EVE_TIMESERIES": "20110607_EVE_L0CS_DIODES_1m.txt",
    "LYRA_LEVEL3_TIMESERIES": "lyra_20110607-000000_lev3_std.fits",
    "GOES_XRS_TIMESERIES": "go1520110607.fits",
    "GBM_TIMESERIES": "glg_cspec_n5_110607_v00.pha",
    "RHESSI_TIMESERIES": "hsi_obssumm_20110607_025.fits",
    "NORH_TIMESERIES": "tca110607.fits",
    "LOFAR_IMAGE": "LOFAR_70MHZ_20190409_131136.fits",
}

# Reverse the dict because we want to use it backwards, but it is nicer to
# write the other way around
_SAMPLE_FILES = {v: k for k, v in _SAMPLE_DATA.items()}


def _download_sample_data(base_url, sample_files, overwrite):
    """
    Downloads a list of files.

    Parameters
    ----------
    base_url : str
        Base URL for each file.
    sample_files : list of tuples
        List of tuples that are (URL_NAME, SAVE_NAME).
    overwrite : bool
        Will overwrite a file on disk if True.

    Returns
    -------
    `parfive.Results`
        Download results. Will behave like a list of files.
    """
    dl = Downloader(overwrite=overwrite, progress=False)
    for url_file_name, fname in sample_files:
        url = urljoin(base_url, url_file_name)
        dl.enqueue_file(url, filename=fname)
    results = dl.download()
    return results


def _retry_sample_data(results):
    # In case we have a broken file on disk, overwrite it.
    dl = Downloader(overwrite=True, progress=False)
    for err in results.errors:
        file_name = err.filepath_partial().name
        log.debug(
            f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        # Update the url to a mirror and requeue the file.
        new_url = urljoin(_BASE_URLS[1], file_name)
        log.debug(f"Attempting redownload of {_SAMPLE_FILES[file_name]} using {new_url}")
        dl.enqueue_file(new_url, filename=err.filepath_partial)
    extra_results = dl.download()
    for err in extra_results.errors:
        file_name = err.filepath_partial().name
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}"
                  )
        log.error(
            f"Failed to download {_SAMPLE_FILES[file_name]} from all mirrors,"
            "the file will not be available."
        )
    return results + extra_results


def download_sample_data(overwrite=False):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    overwrite: `bool`
        Overwrite existing sample data.
    """
    # Workaround for tox only. This is not supported as a user option
    sampledata_dir = os.environ.get("SUNPY_SAMPLEDIR", False)
    if sampledata_dir:
        sampledata_dir = Path(sampledata_dir).expanduser().resolve()
        _is_writable_dir(sampledata_dir)
    else:
        # Creating the directory for sample files to be downloaded
        sampledata_dir = Path(get_and_create_sample_dir())
    already_downloaded = []
    to_download = []
    for url_file_name in _SAMPLE_FILES.keys():
        fname = sampledata_dir/url_file_name
        # We want to avoid calling download if we already have all the files.
        if fname.exists() and not overwrite:
            already_downloaded.append(fname)
        else:
            # URL and Filename pairs
            to_download.append((url_file_name, fname))
    if to_download:
        results = _download_sample_data(_BASE_URLS[0], to_download, overwrite=overwrite)
    else:
        return already_downloaded
    # Something went wrong.
    if results.errors:
        results = _retry_sample_data(results)
    return results + already_downloaded
