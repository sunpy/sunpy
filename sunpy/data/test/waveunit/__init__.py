from pathlib import Path

from sunpy.data.test import rootdir as testrootdir

waveunitdir = str(Path(testrootdir).joinpath('waveunit'))
MEDN_IMAGE = str(Path(waveunitdir).joinpath('medn_halph_fl_20050501_074655.fts'))
MQ_IMAGE = str(Path(waveunitdir).joinpath('mq130812.084253.fits'))
NA_IMAGE = str(Path(waveunitdir).joinpath('na120701.091058.fits'))
SVSM_IMAGE = str(Path(waveunitdir).joinpath('svsm_e3100_S2_20110625_1856.fts'))
