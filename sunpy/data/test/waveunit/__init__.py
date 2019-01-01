import pathlib

from sunpy.data.test import rootdir as testrootdir

waveunitdir = str(pathlib.Path.home().joinpath(testrootdir, 'waveunit'))
MEDN_IMAGE = str(pathlib.Path.home().joinpath(waveunitdir, 'medn_halph_fl_20050501_074655.fts'))
MQ_IMAGE = str(pathlib.Path.home().joinpath(waveunitdir, 'mq130812.084253.fits'))
NA_IMAGE = str(pathlib.Path.home().joinpath(waveunitdir, 'na120701.091058.fits'))
SVSM_IMAGE = str(pathlib.Path.home().joinpath(waveunitdir, 'svsm_e3100_S2_20110625_1856.fts'))
