

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
    "AIA_171_ROLL_IMAGE": "aiacalibim5.fits",
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
    "SRS_TABLE": "20110607SRS.txt",
    "AIA_193_JUN2012": "AIA20120601_000007_0193_lowres.fits",
    "STEREO_A_195_JUN2012": "20120601_000530_n4eua.fits",
    "STEREO_B_195_JUN2012": "20120601_000530_n4eub.fits",
}
