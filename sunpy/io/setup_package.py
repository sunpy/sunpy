from __future__ import absolute_import

import os
import sys
import platform

from distutils.core import Extension
from glob import glob

from astropy_helpers import setup_helpers


def get_extensions():

    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'].append('numpy')
    cfg['sources'].extend(glob(os.path.join(os.path.dirname(__file__), 'src', 'ana', '*.c')))
    if platform.system() == 'Windows':
        cfg['include_dirs'].append(os.path.join(os.path.dirname(__file__), "./msinttypes"))
        cfg['extra_compile_args'].extend(['/D', '"WIN32"',
                                          '/D', '"_WINDOWS"',
                                          '/D', '"_MBCS"',
                                          '/D', '"_USRDLL"',
                                          '/D', '"_CRT_SECURE_NO_DEPRECATE"'])
    else:
        cfg['extra_compile_args'].extend(['-std=c99', '-O3'])
        # Squash some warnings
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
