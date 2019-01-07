# Licensed under a 3-clause BSD style license - see LICENSE.rst
from pathlib import Path


def get_package_data():
    paths = [str(Pathlib.home().joinpath('js', '*.js'))]
    return {'sunpy.extern': paths}
