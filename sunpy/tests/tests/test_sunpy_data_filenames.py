# This test lives here because we have to exclude sunpy/data from pytest.
import os
from pathlib import Path

import sunpy.data.test


def mockreturn(path):
    paths = [
        (os.path.join('test', 'data', ''), (), ('code.py', 'test_file', 'code.pyc', '__init__.py'))
    ]
    return paths


def test_get_test_data_filenames(monkeypatch):
    monkeypatch.setattr(os, 'walk', mockreturn)
    monkeypatch.setattr(os.path, 'isfile', mockreturn)
    output = sunpy.data.test.get_test_data_filenames()
    assert isinstance(output, list)
    # Only the test file and not the py/pyc files should be in the return.
    assert output == [Path('test') / 'data' / 'test_file']
