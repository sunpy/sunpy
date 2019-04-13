"""
The following files are available in this module:
"""
import sys
from pathlib import Path

from ._sample import _sample_files, download_sample_data

files = download_sample_data()

file_list = []
file_dict = {}
for f in files:
    name = Path(f).name
    _key = _sample_files.get(name, None)
    if not _key:
        continue

    setattr(sys.modules[__name__], _key, str(f))
    file_list.append(f)
    file_dict.update({_key: f})
    __doc__ += '* ``{}``\n'.format(_key)

__all__ = list(_sample_files.values()) + ['file_dict', 'file_list']
