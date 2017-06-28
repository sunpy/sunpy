from __future__ import absolute_import

import os
import sys
import platform

from distutils.core import Extension
from glob import glob

from astropy_helpers import setup_helpers


def get_extensions():
    # 'numpy' will be replaced with the proper path to the numpy includes
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'].extend(glob(os.path.join(os.path.dirname(__file__), 'src', 'ana', '*.c')))
    cfg['extra_compile_args'].extend(['-std=c99', '-O3'])

    e = Extension('sunpy.io._pyana', **cfg)
    return [e]

def requires_2to3():
    return False
