import os
import pytest
import sunpy.data.test


@pytest.fixture
def mockreturn(path):
    paths = [
        (os.path.join('test', 'data', ''), (), ('code.py', 'test_file', 'code.pyc', '__init__.py'))
    ]
    return paths


def test_test_data_filenames(monkeypatch):
    monkeypatch.setattr(os, 'walk', mockreturn)
    monkeypatch.setattr(os.path, 'isfile', mockreturn)
    output = sunpy.data.test.test_data_filenames()
    assert isinstance(output, list)
    # Only the test file and not the py/pyc files should be in the return.
    assert output == [os.path.join('test', 'data', '', 'test_file')]
