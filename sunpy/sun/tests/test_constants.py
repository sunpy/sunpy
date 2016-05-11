from __future__ import absolute_import
from sunpy.sun import constants as con
from astropy.constants import Constant
import pytest
from sunpy.extern.six import iteritems


@pytest.mark.parametrize('this_constant', [value for key, value in iteritems(con.constants)])
def test_all_constants_are_constants(this_constant):
    """Test that each member of the constants dict is an astropy Constant"""
    assert type(this_constant) is Constant


@pytest.mark.parametrize('this_key', [key for key, value in iteritems(con.constants)])
def test_get_function(this_key):
    """Test that the get function works for all the keys"""
    assert type(con.get(this_key)) is Constant


@pytest.mark.parametrize('this_key', [key for key, value in iteritems(con.constants)])
def test_find_function(this_key):
    """Test that the find function works for all the keys"""
    assert len(con.find(this_key)) >= 1


@pytest.mark.parametrize('this_key', [key for key, value in iteritems(con.constants)])
def test_find_function(this_key):
    """Test that the find function works for all the keys"""
    assert len(con.find(this_key)) >= 1


@pytest.mark.parametrize("test_input", ['boo', 'crab', 'foo'])
def test_find_function(test_input):
    """Test that the find function fails as expected"""
    assert len(con.find(test_input)) == 0
