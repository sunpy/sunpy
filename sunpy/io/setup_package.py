import os
from glob import glob
from collections import defaultdict

import numpy
from extension_helpers import get_compiler
from setuptools import Extension


def get_extensions():

    if get_compiler() == 'msvc' or os.environ.get("SUNPY_NO_BUILD_ANA_EXTENSION", None):
        return list()
    else:
        cfg = defaultdict(list)
        cfg['include_dirs'].append(numpy.get_include())
        cfg['sources'].extend(sorted(glob(
            os.path.join(os.path.dirname(__file__), 'src', 'ana', '*.c'))))
        cfg['extra_compile_args'].extend(['-std=c99', '-O3'])
        # Squash all warnings
        cfg['extra_compile_args'].extend(['-w'])
        e = Extension('sunpy.io._pyana', **cfg)
        return [e]
