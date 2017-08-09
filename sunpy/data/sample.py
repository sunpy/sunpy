# -*- coding: utf-8 -*-
"""
SunPy sample data files

The following files are available in this submodule:

"""
from __future__ import absolute_import

import sys
from ._sample import _base_urls, _sample_files, get_sample_file

file_list = []
file_dict = {}
for _key in _sample_files:
    f = get_sample_file(_sample_files[_key], _base_urls)
    setattr(sys.modules[__name__], _key, f)
    file_list.append(f)
    file_dict.update({_key: f})
    __doc__ += '* ``{}``\n'.format(_key)

__all__ = list(_sample_files.keys()) + ['file_dict', 'file_list']
