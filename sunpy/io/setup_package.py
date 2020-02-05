import os
from glob import glob
from collections import defaultdict
from distutils.core import Extension

import numpy
from extension_helpers import get_compiler


def get_extensions():

    if get_compiler() == 'msvc':
        return list()
    else:
        cfg = defaultdict(list)
        cfg['include_dirs'].append(numpy.get_include())
        cfg['sources'].extend(sorted(glob(
            os.path.join(os.path.dirname(__file__), 'src', 'ana', '*.c'))))
        cfg['extra_compile_args'].extend(['-std=c99', '-O3'])
        # Squash some warnings
        cfg['extra_compile_args'].extend(['-Wno-unused-but-set-variable',
                                          '-Wno-unused-variable',
                                          '-Wno-unused-result',
                                          '-Wno-sign-compare'])

        e = Extension('sunpy.io._pyana', **cfg)
        return [e]
