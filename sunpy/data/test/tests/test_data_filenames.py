import pytest
import os
import sunpy.data.test


@pytest.fixture
def mockreturn(path):
    return_value = [
        ('/universe', ('solar_system',), ()),
        ('/universe/solar_system', (), ('sun.py', 'earth', 'moon.pyc')),
    ]
    return return_value


def test_test_data_filenames(monkeypatch):
    monkeypatch.setattr(os, 'walk', mockreturn)
    monkeypatch.setattr(os.path, 'isfile', mockreturn)
    output = sunpy.data.test.test_data_filenames()
    all_exist = True
    for file in output:
        if not os.path.isfile(file):
            all_exist = False
            break
    new_op = []
    for path in output:
        new_op.append(path.replace('\\', '/'))  # Handle in case of system = win32
    assert isinstance(new_op, list) and new_op == ['/universe/solar_system/earth'] and all_exist
