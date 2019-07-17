import warnings
from pathlib import Path
from collections import namedtuple
from urllib.parse import urljoin

import parfive

from sunpy.util.config import get_and_create_sample_dir
from sunpy.util.exceptions import SunpyUserWarning

_base_urls = (
    'http://data.sunpy.org/sunpy/v1/',
    'https://github.com/sunpy/sample-data/raw/master/sunpy/v1/',
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
    # Not in the sample-data repo
    # "RHESSI_EVENT_LIST": "hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits",
    "SWAP_LEVEL1_IMAGE": "swap_lv1_20110607_063329.fits",
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

# Reverse the dict because we want to use it backwards, but it is nicer to
# write the other way around
_sample_files = {v: k for k, v in _sample_files.items()}

_error = namedtuple("error", ("filepath_partial", "url", "response"))


def download_sample_data(overwrite=False):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    overwrite: `bool`
        Overwrite existing sample data.
    """
    # Creating the directory for sample files to be downloaded
    sampledata_dir = Path(get_and_create_sample_dir())

    dl = parfive.Downloader(overwrite=overwrite)

    first_url = _base_urls[0]

    already_downloaded = []
    for file_name in _sample_files.keys():
        url = urljoin(first_url, file_name)
        fname = sampledata_dir/file_name
        # We have to avoid calling download if we already have all the files.
        if fname.exists() and not overwrite:
            already_downloaded.append(fname)
        else:
            dl.enqueue_file(url, filename=sampledata_dir/file_name)

    if dl.queued_downloads:
        results = dl.download()
    else:
        return already_downloaded

    if not results.errors:
        return results

    for retry_url in _base_urls[1:]:
        for i, err in enumerate(results.errors):
            file_name = Path(err.url).name
            # Overwrite the parfive error to change the url to a mirror
            new_url = urljoin(retry_url, file_name)
            results._errors[i] = _error(err.filepath_partial,
                                        new_url,
                                        err.exception)

        results = dl.retry(results)

        if not results.errors:
            return results

    for err in results.errors:
        file_name = Path(err.url).name
        warnings.warn(f"File {file_name} not found.", SunpyUserWarning)

    return results
