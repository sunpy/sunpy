"""
This module provides the following sample data files.  These files are
downloaded when this module is imported for the first time.

.. list-table::
   :widths: auto
   :header-rows: 1

   * - Variable name
     - Name of downloaded file
"""
import sys
from pathlib import Path

from ._sample import _SAMPLE_FILES, download_sample_data

files = download_sample_data()

file_dict = {}
for f in files:
    name = Path(f).name
    _key = _SAMPLE_FILES.get(name, None)
    if _key:
        setattr(sys.modules[__name__], _key, str(f))
        file_dict.update({_key: f})

# Sort the entries in the dictionary
file_dict = dict(sorted(file_dict.items()))

file_list = file_dict.values()

for keyname, filename in file_dict.items():
    __doc__ += f'   * - ``{keyname}``\n     - {Path(filename).name}\n'

__all__ = list(_SAMPLE_FILES.values()) + ['file_dict', 'file_list']
