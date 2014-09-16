# -*- coding: utf-8 -*-
"""SunPy sample data files"""

import sys
from _sample import sample_files
from os.path import isfile, join

from sunpy import config
sampledata_dir = config.get("downloads", "sample_dir")

for key in sample_files:
    print(join(sampledata_dir, sample_files[key]))
    if isfile(join(sampledata_dir, sample_files[key])):
        setattr(sys.modules[__name__], key, join(sampledata_dir, sample_files[key]))
    else:
        raise ImportError("Sample data file(s) missing. Use download() to get them.")