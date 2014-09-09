from __future__ import absolute_import

from sunpy.data import sample

for key in sample.sample_files:
    setattr(sample, key, sample.sample_files[key])