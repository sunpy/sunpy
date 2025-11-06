import os
from pathlib import Path
from collections import defaultdict

import numpy
from extension_helpers import get_compiler
from setuptools import Extension


def get_extensions():
    if os.environ.get("SUNPY_NO_BUILD_ANA_EXTENSION"):
        return []
    cfg = defaultdict(list)
    cfg["include_dirs"].append(numpy.get_include())
    cfg["sources"].extend(list((Path(__file__).parent / "src" / "ana").glob("*.c")))
    if get_compiler() == 'msvc':
        cfg["extra_compile_args"].extend([
            "/O2",
            "/W3",
            "/utf-8",
        ])
        # quiet MSVC's CRT nags
        cfg["define_macros"].extend([
            ("_CRT_SECURE_NO_WARNINGS", None),
        ])
    else:
        cfg["extra_compile_args"].extend([
            "-std=c99",
            "-O3",
            "-w",
        ])

    return [Extension("sunpy.io._pyana", **cfg)]
