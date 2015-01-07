"""This module tests the functions implemented in sunpy.util.util."""

from __future__ import absolute_import

from sunpy.util import util
import numpy as np
import warnings

def test_to_signed():
    """
    This should return a signed type that can hold uint32.
    """
    assert util.to_signed(np.dtype('uint32')) == np.dtype('int64')

def test_goes_flare_class():
    """
    This should convert a list of GOES classes into a list of numbers.
    """
    lst = ['A1.0', 'M1.0', 'C1.0', 'B1.0', 'X1.0']
    assert util.goes_flare_class(lst) == \
        [1.0e-08, 1.0e-05, 1.0e-06, 1.0e-07, 1.0e-04]

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
    lst = [['n', 'sqrt(n)', 'n^2'], \
           ['1', '1', '1'], \
           ['4', '2', '16'], \
           ['3', '1.732', '9']]
    expected = ('n|sqrt(n)|n^2\n'
               '1|1      |1  \n'
               '4|2      |16 \n'
               '3|1.732  |9  ')
    assert util.print_table(lst, colsep='|') == expected

def test_findpeaks():
    """
    This should return the indices of the local maxima of numpy array
    data (relative to index 1).
    """
    data = np.array([1.0, 3.5, 3.0, 4.0, -9.0, 0.0, 0.5, 0.3, 9.5])
    assert np.array_equal(util.findpeaks(data), np.array([0, 2, 5]))

def test_polyfun_at():
    """
    This should evaluate the polynomial x^3 + 5x^2 - 6x + 3 at x = 5.
    """
    coeff = [1, 5, -6, 3]
    assert util.polyfun_at(coeff, 5) == 223

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
    assert list(util.find_next(list1, list2, None)) == [(1, 2), (2, 3), (3, 5),
                        (4, 5), (5, 9), (6, 10), (7, 15), (8, None), (9, None)]

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

def test_replacement_filename():
    """
    This should return a replacement path for the current file.
    """
    assert util.replacement_filename(__file__).endswith('test_util.0.py')

def test_expand_list():
    """
    This should return an expanded version of list lst.
    """
    lst = [1, 2, 3, [4, 5, 6], 7, (8, 9)]
    assert util.expand_list(lst) == [1, 2, 3, 4, 5, 6, 7, 8, 9]

def test_deprecated():
    """
    This should trigger a deprecation warning.
    """
    depr = util.Deprecated()
    with warnings.catch_warnings(record=True) as current_warnings:
        depr_func = depr(lambda x: x)
        depr_func(1)
        assert len(current_warnings) == 1
