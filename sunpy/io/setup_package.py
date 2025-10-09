import os
import sys
from glob import glob
from collections import defaultdict

import numpy
from setuptools import Extension


def get_extensions():
    if os.environ.get("SUNPY_NO_BUILD_ANA_EXTENSION"):
        return []

    cfg = defaultdict(list)
    cfg["include_dirs"].append(numpy.get_include())
    cfg["sources"].extend(
        sorted(glob(os.path.join(os.path.dirname(__file__), "src", "ana", "*.c")))
    )
    # quiet MSVC's CRT nags; harmless elsewhere if ignored
    cfg["define_macros"].extend([
        ("_CRT_SECURE_NO_WARNINGS", None),
        ("_CRT_NONSTDC_NO_DEPRECATE", None),
    ])

    if sys.platform.startswith("win"):
        cfg["extra_compile_args"].extend([
            "/O2",
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
