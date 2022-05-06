#!/bin/env python

"""
This script is used to update the files in `sunpy/extern` while making new releases.
"""

# TODO:
# [x] Get list of all the packages imported in sunpy/extern
# [] Search for the source distribution of each package in `USER/home/sunpy/extern`
# [] Extract the source distribution to a temporary directory
# [] Update the python files in `sunpy/extern`


import ast
from pathlib import Path

# get parent directory of sunpy
SUNPY_DIR = Path(__file__).parent.parent


def visit_import(node):
    for name in node.names:
        modules.append(name.name.split(".")[0])


def visit_from(node):
    # if node.module is missing it's a "from . import ..." statement
    # if level > 0 it's a "from .submodule import ..." statement
    if node.module is not None and node.level == 0:
        modules.append(node.module.split(".")[0])


# get all the modules imported in sunpy/extern
file_names = ["appdirs", "distro", "inflect", "parse"]

modules = list()

for i in file_names:
    with open(SUNPY_DIR/f"sunpy/extern/{i}.py", "r") as f:
        # print(f"Reading {i}.py")
        tree = ast.parse(f.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                visit_import(node)
            elif isinstance(node, ast.ImportFrom):
                visit_from(node)

# remove python modules from the list
modules = [i for i in modules if i not in ['platform', 'warnings', 're', 'ctypes', 'os',
                                           'shlex', 'functools', 'ast',
                                           'datetime', 'logging', 'sys', 'com', 'decimal', 'argparse', '__future__', 'json', 'typing', 'subprocess', 'array']]
modules = set(modules)
print(modules)

# Output:
# {'win32api', 'win32com', 'winreg', '_winreg'}
