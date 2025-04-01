"""
This module provides the following sample data files. When a sample shortname
is accessed, the corresponding file is downloaded if needed. All files can be
downloaded by calling :func:`~sunpy.data.sample.download_all`.

Summary variables
-----------------
.. list-table::
   :widths: auto

   * - ``file_dict``
     - Dictionary of all sample shortnames and, if downloaded, corresponding
       file locations on disk (otherwise, ``None``)
   * - ``file_list``
     - List of disk locations for sample data files that have been downloaded

Sample shortnames
-----------------
.. list-table::
   :widths: auto
   :header-rows: 1

   * - Sample shortname
     - Name of downloaded file
"""
from ._sample import _SAMPLE_DATA, _get_sample_files

# Add a table row to the module docstring for each sample file
for _keyname, _filename in sorted(_SAMPLE_DATA.items()):
    __doc__ += f'   * - ``{_keyname}``\n     - {_filename}\n'


# file_dict and file_list are not normal variables; see __getattr__() below
__all__ = list(sorted(_SAMPLE_DATA.keys())) + ['download_all', 'file_dict', 'file_list']  # NOQA: F822


# See PEP 562 (https://peps.python.org/pep-0562/) for module-level __dir__()
def __dir__():
    return __all__


# See PEP 562 (https://peps.python.org/pep-0562/) for module-level __getattr__()
def __getattr__(name):
    if name in _SAMPLE_DATA:
        return _get_sample_files([_SAMPLE_DATA[name]])[0]
    elif name == 'file_dict':
        return dict(sorted(zip(_SAMPLE_DATA.keys(),
                               _get_sample_files(_SAMPLE_DATA.values(), no_download=True))))
    elif name == 'file_list':
        return [v for v in __getattr__('file_dict').values() if v]
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def download_all(force_download=False):
    """
    Download all sample data at once that has not already been downloaded.

    Parameters
    ----------
    force_download : `bool`
        If ``True``, files are downloaded even if they already exist. Default is
        ``False``.
    """
    _get_sample_files(_SAMPLE_DATA.values(), force_download=force_download)
