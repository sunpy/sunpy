import pytest

from astropy.constants import Constant
from astropy.table import Table

from sunpy.sun import constants as con


def test_find_all():
    assert isinstance(con.find(), list)
    assert len(con.find()) == 34


def test_print_all():
    table = con.print_all()
    assert isinstance(table, Table)
    assert len(table) == 34


@pytest.mark.parametrize('this_constant', [value for key, value in con.constants.items()])
def test_all_constants_are_constants(this_constant):
    """
    Test that each member of the constants dict is an astropy Constant.
    """
    assert isinstance(this_constant, Constant)


@pytest.mark.parametrize('this_key', [key for key, value in con.constants.items()])
def test_get_function(this_key):
    """
    Test that the get function works for all the keys.
    """
    assert isinstance(con.get(this_key), Constant)


@pytest.mark.parametrize('this_key', [key for key, value in con.constants.items()])
def test_find_function(this_key):
    """
    Test that the find function works for all the keys.
    """
    assert len(con.find(this_key)) >= 1


@pytest.mark.parametrize('this_key', [key for key, value in con.constants.items()])
def test_find_function2(this_key):
    """
    Test that the find function works for all the keys.
    """
    assert len(con.find(this_key)) >= 1


@pytest.mark.parametrize("test_input", ['boo', 'crab', 'foo'])
def test_find_function3(test_input):
    """
    Test that the find function fails as expected.
    """
    assert len(con.find(test_input)) == 0


def test_docstring():
    """
    Test that the docstring RST table has the correct number of constants.
    """
    lines = con.__doc__.split('\n')
    description = 'The following constants are available:'
    assert description in lines

    # After the description line, there are five lines before the actual table
    # data begins (including a newline, RST headers, and column names). Count
    # the number of rows in the table until the RST column footer is reached.
    data_start_idx = lines.index(description) + 5
    data_end_idx = data_start_idx
    for idx, line in enumerate(lines[data_start_idx:], data_start_idx):
        if line.startswith('='):
            data_end_idx = idx
            break

    num_rows = data_end_idx - data_start_idx
    assert num_rows == 34
