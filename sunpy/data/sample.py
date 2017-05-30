# -*- coding: utf-8 -*-
"""SunPy sample data files"""
from __future__ import absolute_import
from sunpy.extern import six

import sys
from ._sample import sample_files as _sample_files
from ._sample import download_sample_data
import os.path

from sunpy import config as _config
_sampledata_dir = _config.get("downloads", "sample_dir")

_base_urls = (
    'http://data.sunpy.org/sample-data/',
    'https://github.com/sunpy/sunpy-sample-data/raw/master/')

for _key in _sample_files:    # remove zip extension if exists
    if _sample_files[_key][-3:] == 'zip':
        f = _sample_files[_key][:-4]
    else:
        f = _sample_files[_key]
    if not os.path.isfile(os.path.join(_sampledata_dir, f)):
        download_sample_data(f, _base_urls)
    setattr(sys.modules[__name__], _key, os.path.join(_sampledata_dir, f))

file_dict = _sample_files
file_list = list(_sample_files.values())

__all__ = list(_sample_files.keys()) + ['file_dict', 'file_list']
