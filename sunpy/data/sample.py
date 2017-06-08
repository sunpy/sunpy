# -*- coding: utf-8 -*-
"""SunPy sample data files"""
from __future__ import absolute_import

import sys
from ._sample import _files as _sample_files
from ._sample import get_sample_file

_base_urls = (
    'http://data.sunpy.org/sample-data/',
    'https://github.com/sunpy/sunpy-sample-data/raw/master/'
)

file_list = []
file_dict = {}
for _key in _sample_files:
    f = get_sample_file(_sample_files[_key], _base_urls)
    setattr(sys.modules[__name__], _key, f)
    file_list.append(f)
    file_dict.update({_key: f})

__all__ = list(_sample_files.keys()) + ['file_dict', 'file_list']
