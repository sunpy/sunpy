from __future__ import absolute_import

import os
import platform

from distutils.core import Extension
from glob import glob

from astropy_helpers import setup_helpers


def get_extensions():
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'].extend(glob(os.path.join(os.path.dirname(__file__), 'src', 'ana', '*.c')))
    if platform.system() == 'Windows' and int(platform.python_version()[0]) > 2:
        cfg['include_dirs'].append(os.path.join(os.path.dirname(__file__), "msinttypes"))
    elif platform.system() == 'Windows' and int(platform.python_version()[0]) == 2:
        return list()
    else:
        cfg['extra_compile_args'].extend(['-std=c99', '-O3'])
        # Squash ALL (probably) warnings
        cfg['extra_compile_args'].extend(['-Wno-declaration-after-statement',
                                          '-Wno-unused-variable', '-Wno-parentheses',
                                          '-Wno-uninitialized', '-Wno-format',
                                          '-Wno-strict-prototypes', '-Wno-unused', '-Wno-comments',
                                          '-Wno-switch', '-Wno-strict-aliasing', '-Wno-return-type',
                                          '-Wno-address'])
    e = Extension('sunpy.io._pyana', **cfg)
    return [e]

def requires_2to3():
    return False
