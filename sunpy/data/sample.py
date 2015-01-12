# -*- coding: utf-8 -*-
"""SunPy sample data files"""

import sys as _sys
from _sample import sample_files as _sample_files
from os.path import isfile as _isfile
from os.path import join as _join

from sunpy import config as _config
_sampledata_dir = _config.get("downloads", "sample_dir")

for key in _sample_files:
    if _isfile(_join(_sampledata_dir, _sample_files[key])):
        setattr(_sys.modules[__name__], key, _join(_sampledata_dir, _sample_files[key]))
    else:
        raise ImportError("Sample data file(s) missing. Use sunpy.data.download_sample_data() to get them.")