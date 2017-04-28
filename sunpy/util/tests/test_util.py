"""
This module tests the functions implemented in sunpy.util.util.
"""

from __future__ import absolute_import, division, print_function
import numpy as np

import pytest
from pytest_mock import mocker

from sunpy.util import util


def test_to_signed():
    """
    This should return a signed type that can hold uint32 and ensure that
    an exception is raised when attempting to convert an unsigned 64 bit integer
    to an integer
    """
    assert util.to_signed(np.dtype('uint32')) == np.dtype('int64')

    with pytest.raises(ValueError):
        util.to_signed(np.dtype('uint64')) == np.dtype('int64')


def test_unique():
    """
    This should add the unique values of itr to unique_list.
    """
    itr = [6, 1, 2, 1, 7, 41.2, '41.2', 1, '41.2']
    unique_list = []
    for elem in util.unique(itr):
        unique_list.append(elem)
    assert unique_list == [6, 1, 2, 7, 41.2, '41.2']


def test_unique_key():
    """
    This should add each element of itr to unique_list if no preceding
    element is congruent to it in mod 10.
    """
    itr = [7, 3, 17, 104, 6, 1006, 117, 14, 10]
    unique_list = []
    for elem in util.unique(itr, lambda x: x % 10):
        unique_list.append(elem)
    assert unique_list == [7, 3, 104, 6, 10]


def test_print_table():
    """
    This should return a string representation of lst with table elements
    left-justified and with columns separated by dashes.
    """
    lst = [['n', 'sqrt(n)', 'n^2'],
           ['1', '1', '1'],
           ['4', '2', '16'],
           ['3', '1.732', '9']]
    expected = ('n|sqrt(n)|n^2\n'
                '1|1      |1  \n'
                '4|2      |16 \n'
                '3|1.732  |9  ')
    assert util.print_table(lst, colsep='|') == expected


def test_minimal_pairs():
    """
    This should return the pairs of elements from list1 and list2 with
    minimal difference between their values.
    """
    list1 = [0, 5, 10, 15, 20, 25]
    list2 = [3, 12, 19, 21, 26, 29]
    assert list(util.minimal_pairs(list1, list2)) == [(1, 0, 2), (2, 1, 2),
                                                      (4, 2, 1), (5, 4, 1)]


def test_find_next():
    """
    This should return a generator yielding the nearest larger element in
    list2 for each element in list1 (or None if none exists after the
    previous element yielded from list2).
    """
    list1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    list2 = [0, 2, 3, 5, 0, 0, 5, 9, 10, 15]
    assert list(util.find_next(list1, list2, None)) == [(1, 2), (2, 3), (3, 5), (4, 5), (5, 9),
                                                        (6, 10), (7, 15), (8, None), (9, None)]


def test_common_base():
    """
    This should return the base class common to each object in objs.
    """
    class TestA(object):
        """Base test class."""
        pass

    class TestB(TestA):
        """First inherited class."""
        pass

    class TestC(TestA):
        """Second inherited class."""
        pass
    inst_b = TestB()
    inst_c = TestC()
    objs = [inst_b, inst_c]
    assert util.common_base(objs) == TestA


def test_merge():
    """
    This should return a sorted (from greatest to least) merged list
    from list1 and list2.
    """
    list1 = [13, 11, 9, 7, 5, 3, 1]
    list2 = [14, 12, 10, 8, 6, 4, 2]
    result = list(util.merge([list1, list2]))
    assert result[::-1] == sorted(result)

    assert list(util.merge([[], [1], []])) == [1]


def test_replacement_filename():
    """
    This should return a replacement path for the current file.
    """
    assert util.replacement_filename(__file__).endswith('test_util.0.py')


def test_replacement_filename_path_not_exists(mocker):
    """
    If a candidate path does not exist, then just return it as it is OK to use
    """
    path_not_exists = '/tmp'
    mocker.patch('os.path.exists', return_value=False)

    assert util.replacement_filename(path_not_exists) == path_not_exists


def test_expand_list():
    """
    This should return an expanded version of list lst.
    """
    lst = [1, 2, 3, [4, 5, 6], 7, (8, 9)]
    assert util.expand_list(lst) == [1, 2, 3, 4, 5, 6, 7, 8, 9]


def test_expand_list_generator():

    lst = ['a', 'b', [], (['c', 'd']), tuple(), ['e']]
    assert list(util.expand_list_generator(lst)) == ['a', 'b', 'c', 'd', 'e']
