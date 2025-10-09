import os
from glob import glob
from collections import defaultdict

import numpy
from extension_helpers import get_compiler
from setuptools import Extension


def get_extensions():
    if os.environ.get("SUNPY_NO_BUILD_ANA_EXTENSION"):
        return []
    cfg = defaultdict(list)
    cfg["include_dirs"].append(numpy.get_include())
    cfg["sources"].extend(
        sorted(glob(os.path.join(os.path.dirname(__file__), "src", "ana", "*.c")))
    )

    if get_compiler() == 'msvc':
        cfg["extra_compile_args"].extend([
            "/O3",
            "/W3",
            "/utf-8",
        ])
    else:
        cfg["extra_compile_args"].extend([
            "-std=c99",
            "-O3",
            "-w",
        ])

    return [Extension("sunpy.io._pyana", **cfg)]
