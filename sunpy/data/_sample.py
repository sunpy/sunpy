import os
from pathlib import Path
from urllib.parse import urljoin

from sunpy import log
from sunpy.util.config import _is_writable_dir, get_and_create_sample_dir
from sunpy.util.parfive_helpers import Downloader

_BASE_URLS = (
    'https://github.com/sunpy/data/raw/main/sunpy/v1/',
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
    "AIA_094_IMAGE": "AIA20110607_063305_0094_lowres.fits",
    "AIA_131_IMAGE": "AIA20110607_063301_0131_lowres.fits",
    "AIA_1600_IMAGE": "AIA20110607_063305_1600_lowres.fits",
    "AIA_1600_VENUS_IMAGE": "aia_lev1_1600a_2012_06_06t04_07_29_12z_image_lev1_lowres.fits",
    "AIA_171_IMAGE": "AIA20110607_063302_0171_lowres.fits",
    "AIA_171_ROLL_IMAGE": "aiacalibim5.fits",
    "AIA_193_CUTOUT01_IMAGE": "AIA20110607_063307_0193_cutout.fits",
    "AIA_193_CUTOUT02_IMAGE": "AIA20110607_063931_0193_cutout.fits",
    "AIA_193_CUTOUT03_IMAGE": "AIA20110607_064555_0193_cutout.fits",
    "AIA_193_CUTOUT04_IMAGE": "AIA20110607_065219_0193_cutout.fits",
    "AIA_193_CUTOUT05_IMAGE": "AIA20110607_065843_0193_cutout.fits",
    "AIA_193_IMAGE": "AIA20110607_063307_0193_lowres.fits",
    "AIA_193_JUN2012": "AIA20120601_000007_0193_lowres.fits",
    "AIA_211_IMAGE": "AIA20110607_063302_0211_lowres.fits",
    "AIA_304_IMAGE": "AIA20110607_063334_0304_lowres.fits",
    "AIA_335_IMAGE": "AIA20110607_063303_0335_lowres.fits",
    "CALLISTO_SPECTRUM": "BIR_20110607_062400_10.fit",
    "EIT_195_IMAGE": "eit_l1_20110607_203753.fits",
    "EVE_TIMESERIES": "20110607_EVE_L0CS_DIODES_1m.txt",
    "GBM_TIMESERIES": "glg_cspec_n5_110607_v00.pha",
    "GOES_XRS_TIMESERIES": "go1520110607.fits",
    "HMI_LOS_IMAGE": "HMI20110607_063211_los_lowres.fits",
    "LOFAR_IMAGE": "LOFAR_70MHZ_20190409_131136.fits",
    "LYRA_LEVEL3_TIMESERIES": "lyra_20110607-000000_lev3_std.fits",
    "NORH_TIMESERIES": "tca110607.fits",
    "RHESSI_IMAGE": "hsi_image_20110607_063300.fits",
    "RHESSI_TIMESERIES": "hsi_obssumm_20110607_025.fits",
    "SRS_TABLE": "20110607SRS.txt",
    "STEREO_A_195_JUN2012": "20120601_000530_n4eua.fits",
    "STEREO_B_195_JUN2012": "20120601_000530_n4eub.fits",
    "SWAP_LEVEL1_IMAGE": "swap_lv1_20110607_063329.fits",
    "AIA_STEREOSCOPIC_IMAGE": "aia_lev1_171a_2023_07_06t00_05_33_35z_image_lev1.fits",
    "EUVI_STEREOSCOPIC_IMAGE": "20230706_000525_n4eua.fts"
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
    dl = Downloader(overwrite=overwrite, progress=True)

    for url_file_name, fname in sample_files:
        url = urljoin(base_url, url_file_name)
        dl.enqueue_file(url, filename=fname)

    results = dl.download()
    return results


def _retry_sample_data(results, new_url_base):
    # In case we have a broken file on disk, overwrite it.
    dl = Downloader(overwrite=True, progress=True)

    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(
            f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        # Update the url to a mirror and requeue the file.
        new_url = urljoin(new_url_base, file_name)
        log.debug(f"Attempting redownload of {_SAMPLE_FILES[file_name]} using {new_url}")
        dl.enqueue_file(new_url, filename=err.filepath_partial)

    extra_results = dl.download()

    # Make a new results object which contains all the successful downloads
    # from the previous results object and this retry, and all the errors from
    # this retry.
    new_results = results + extra_results
    new_results._errors = extra_results._errors
    return new_results


def _handle_final_errors(results):
    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}"
                  )
        log.error(
            f"Failed to download {_SAMPLE_FILES[file_name]} from all mirrors,"
            "the file will not be available."
        )


def _get_sampledata_dir():
    # Workaround for tox only. This is not supported as a user option
    sampledata_dir = os.environ.get("SUNPY_SAMPLEDIR", False)
    if sampledata_dir:
        sampledata_dir = Path(sampledata_dir).expanduser().resolve()
        _is_writable_dir(sampledata_dir)
    else:
        # Creating the directory for sample files to be downloaded
        sampledata_dir = Path(get_and_create_sample_dir())
    return sampledata_dir


def _get_sample_files(filename_list, no_download=False, force_download=False):
    """
    Returns a list of disk locations corresponding to a list of filenames for
    sample data, downloading the sample data files as necessary.

    Parameters
    ----------
    filename_list : `list` of `str`
        List of filenames for sample data
    no_download : `bool`
        If ``True``, do not download any files, even if they are not present.
        Default is ``False``.
    force_download : `bool`
        If ``True``, download all files, and overwrite any existing ones.
        Default is ``False``.

    Returns
    -------
    `list` of `pathlib.Path`
        List of disk locations corresponding to the list of filenames. An entry
        will be ``None`` if ``no_download == True`` and the file is not present.

    Raises
    ------
    RuntimeError
        Raised if any of the files cannot be downloaded from any of the mirrors.
    """
    sampledata_dir = _get_sampledata_dir()

    fullpaths = [sampledata_dir/fn for fn in filename_list]

    if no_download:
        fullpaths = [fp if fp.exists() else None for fp in fullpaths]
    else:
        to_download = zip(filename_list, fullpaths)
        if not force_download:
            to_download = [(fn, fp) for fn, fp in to_download if not fp.exists()]

        if to_download:
            results = _download_sample_data(_BASE_URLS[0], to_download, overwrite=force_download)

            # Try the other mirrors for any download errors
            if results.errors:
                for next_url in _BASE_URLS[1:]:
                    results = _retry_sample_data(results, next_url)
                    if not results.errors:
                        break
                else:
                    _handle_final_errors(results)
                    raise RuntimeError

    return fullpaths
