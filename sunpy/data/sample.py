# -*- coding: utf-8 -*-
"""SunPy sample data files"""

import sys
from _sample import sample_files as _sample_files
import os.path

from sunpy import config as _config
_sampledata_dir = _config.get("downloads", "sample_dir")

for _key in _sample_files:
    if os.path.isfile(os.path.join(_sampledata_dir, _sample_files[_key])):
        setattr(sys.modules[__name__], _key, os.path.join(_sampledata_dir, _sample_files[_key]))
    else:
        raise ImportError("Sample data file(s) missing. Use sunpy.data.download_sample_data() to get them.")

file_dict = _sample_files
file_list = _sample_files.values()

__all__ = _sample_files.keys() + ['file_dict', 'file_list']
