"""
``sunpy``
=========

An open-source Python library for Solar Physics data analysis.

* Homepage: https://sunpy.org
* Documentation: https://docs.sunpy.org/en/stable/
"""
from sunpy.tests.self_test import self_test
from sunpy.util import system_info
from sunpy.util.config import load_config, print_config
from sunpy.util.logger import _init_log
from .version import version as __version__


def _get_bibtex():
    import os
    import textwrap

    # Set the bibtex entry to the article referenced in CITATION.rst
    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION.rst')

    # Explicitly specify UTF-8 encoding in case the system's default encoding is problematic
    with open(citation_file, encoding='utf-8') as citation:
        # Extract the first bibtex block:
        ref = citation.read().partition(".. code:: bibtex\n\n")[2]
        lines = ref.split("\n")
        # Only read the lines which are indented
        lines = lines[:[line.startswith("    ") for line in lines].index(False)]
        ref = textwrap.dedent('\n'.join(lines))
    return ref


__citation__ = __bibtex__ = _get_bibtex()
config = load_config()
log = _init_log(config=config)

__all__ = ['config', 'self_test', 'system_info', 'print_config', 'log', '__version__', '__citation__', '__bibtex__']
