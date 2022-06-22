"""
This module provides the following sample data files.  These files are
downloaded when this module is imported for the first time.  See
:ref:`sphx_glr_generated_gallery_acquiring_data_2011_06_07_sampledata_overview.py`
for plots of some of these files.

.. list-table::
   :widths: auto
   :header-rows: 1

   * - Variable name
     - Name of downloaded file
"""
import os
from pathlib import Path

import pooch

from sunpy.util.config import _is_writable_dir, get_and_create_sample_dir
from ._sample import _SAMPLE_DATA


def get_sample_data_dir() -> Path:
    # Workaround for tox only. This is not supported as a user option
    sampledata_dir = os.environ.get("SUNPY_SAMPLEDIR", False)
    if sampledata_dir:
        sampledata_dir = Path(sampledata_dir).expanduser().resolve()
        _is_writable_dir(sampledata_dir)
    else:
        # Creating the directory for sample files to be downloaded
        sampledata_dir = Path(get_and_create_sample_dir())
    return sampledata_dir


AVAILABLE_DATA = list(_SAMPLE_DATA.keys())
FILES = list(_SAMPLE_DATA.values())
REGISTRY = {fname: None for fname in FILES}
version = 'v1'
POOCH = pooch.create(
    # Use the default cache folder for the operating system
    path=get_sample_data_dir(),
    # The remote data is on Github
    base_url="https://github.com/sunpy/sample-data/raw/master/sunpy/{version}/",
    version='v1',
    registry=REGISTRY,
)


def __getattr__(value):
    if value in ['file_dict', 'file_list']:
        data_dir = get_sample_data_dir()
        file_dict = {
            keyname: file_path for keyname, filename in _SAMPLE_DATA.items()
            if (file_path := data_dir / filename).exists()
        }
        file_dict = dict(sorted(file_dict.items()))
        if value == 'file_dict':
            return file_dict
        elif value == 'file_list':
            return file_dict.values()
    if value not in AVAILABLE_DATA:
        raise AttributeError(f'{value} not in sample data set:\n{AVAILABLE_DATA}')
    return POOCH.fetch(_SAMPLE_DATA[value])


def download_sample_data(overwrite=False):
    return [POOCH.fetch(fname) for fname in FILES]


for keyname, filename in _SAMPLE_DATA.items():
    __doc__ += f'   * - ``{keyname}``\n     - {filename}\n'

__all__ = AVAILABLE_DATA
