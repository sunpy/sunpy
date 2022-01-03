import os
from glob import glob
from collections import defaultdict

import numpy
from extension_helpers import get_compiler
from setuptools import Extension


def get_extensions():
    exts = []

    if get_compiler() == 'msvc':
        return exts

    try:
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

        exts.append(Extension('sunpy.io._pyana', **cfg))
    except Exception as exc:
        log.info(f"Failed to compile sunpy.io.ana with the following error:\n{exc}")
    finally:
        return exts
