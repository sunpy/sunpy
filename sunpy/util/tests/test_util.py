from inspect import cleandoc

import numpy as np
import pytest

from astropy.io import fits

from sunpy.data.test import get_test_filepath
from sunpy.util import util
from sunpy.util.parfive_helpers import Results


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
    This should add each element of itr to unique_list if no preceding element
    is congruent to it in mod 10.
    """
    itr = [7, 3, 17, 104, 6, 1006, 117, 14, 10]
    unique_list = []
    for elem in util.unique(itr, lambda x: x % 10):
        unique_list.append(elem)
    assert unique_list == [7, 3, 104, 6, 10]


def test_replacement_filename():
    """
    This should return a replacement path for the current file.
    """
    assert util.replacement_filename(__file__).endswith('test_util.0.py')


def test_replacement_filename_path_not_exists(mocker):
    """
    If a candidate path does not exist, then just return it as it is OK to use.
    """
    path_not_exists = '/tmp'
    mocker.patch('os.path.exists', return_value=False)

    assert util.replacement_filename(path_not_exists) == path_not_exists


def test_expand_list():
    """
    This should return an expanded version of list lst.
    """
    lst = [1, 2, 3, [4, 5, 6], 7, (8, 9), ((10, 11), ((12, 13),))]
    assert util.expand_list(lst) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]


@pytest.mark.parametrize(('input_data', 'expected_output'), [
    (['a1234', 'b', [], (['c', 'd']), (), ['e'], b'fghj'], ['a1234', 'b', 'c', 'd', 'e', b'fghj']),
    (iter(['c', 'd']), ['c', 'd']),
    (zip(['1', '2'], ['3', '4']), ['1', '3', '2', '4']),
    (Results(["a", "b", "c"]), ["a", "b", "c"]),
    ([open(get_test_filepath('aia_171_level1.fits'), 'rb')], open(get_test_filepath('aia_171_level1.fits'), 'rb').readlines())
])
def test_expand_list_generator(input_data, expected_output):
    assert list(util.expand_list_generator(input_data)) == expected_output


@pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")
def test_expand_list_generator_map(aia171_test_map):
    before = after = [aia171_test_map.data, aia171_test_map.wcs, fits.getheader(get_test_filepath('aia_171_level1.fits'))]
    assert list(util.expand_list_generator(before)) == after


def test_partial_key_match():
    test_dict = {('a', 'b', 'c'): (1, 2, 3), ('a', 'b', 'd'): (4, 5, 6), ('e', 'f', 'g'): (8, 7, 9)}
    assert list(util.partial_key_match(('a', None, 'c'), test_dict))[
        0] == test_dict[('a', 'b', 'c')]


def test_dict_keys_same():
    dicts = [{'x': 42}, {'x': 23, 'y': 5}]
    assert util.dict_keys_same(dicts) == [{'y': None, 'x': 42}, {'y': 5, 'x': 23}]


def test_get_keywords():

    def f(a, b, c=1, d=2, **e):
        pass

    def g(a, b, *, c=1, h=2, **e):
        pass

    assert util.get_keywords(f) == {'c', 'd'}  # POSITIONAL_OR_KEYWORD
    assert util.get_keywords(g) == {'c', 'h'}  # KEYWORD_ONLY
    assert util.get_keywords([f, g]) == {'c', 'd', 'h'}


def test_get_set_methods():

    class A:

        def _set_test1(self, *args, **kwargs):
            pass

        def set_test2(self, *args, **kwargs):
            pass

        @property
        def set_test3(self):
            pass

        def test4(self, *args, **kwargs):
            pass

        @property
        def test5(self, *args, **kwargs):
            pass

        def set_test6_n(self, *args, **kwargs):
            pass

    assert util.get_set_methods(A()) == {'test2', 'test6_n'}


def test_notes_combined():
    original_documentation = """
    Class Info.

    Notes
    -----
    This is a note.

    References
    ----------
    This is reference.
    """
    extra_note_section = """\nNotes\n-----\nThis should be combined."""
    updated_documentation = util.fix_duplicate_notes(extra_note_section, original_documentation)
    expected_result = """
    Class Info.

    Notes
    -----
    This is a note.

    This should be combined.

    References
    ----------
    This is reference.
    """
    assert updated_documentation == cleandoc(expected_result)


def test_notes_combined_no_references():
    original_documentation = """
    Class Info.

    Notes
    -----
    This is a note.
    """
    extra_note_section = """\nNotes\n-----\nThis should be combined."""
    updated_documentation = util.fix_duplicate_notes(extra_note_section, original_documentation)
    expected_result = """
    Class Info.

    Notes
    -----
    This is a note.

    This should be combined.
    """
    assert updated_documentation == cleandoc(expected_result)


def test_notes_combined_no_existing_notes():
    original_documentation = """
    Class Info.

    References
    ----------
    This is reference.
    """
    extra_note_section= """\nNotes\n-----\nThis should be combined."""
    updated_documentation = util.fix_duplicate_notes(extra_note_section, original_documentation)
    expected_result = """
    Class Info.

    Notes
    -----
    This should be combined.

    References
    ----------
    This is reference.
    """
    assert updated_documentation == cleandoc(expected_result)


def test_notes_combined_no_notes_no_references():
    original_documentation = """
    Class Info.
    """
    extra_note_section= """\nNotes\n-----\nThis should be combined."""
    updated_documentation= util.fix_duplicate_notes(extra_note_section, original_documentation)
    expected_result = """
    Class Info.

    Notes
    -----
    This should be combined.
    """
    assert updated_documentation == cleandoc(expected_result)


def test_grid_perimeter():
    assert np.all(util.grid_perimeter(2, 3) == [[0, 0], [0, 1], [0, 2], [0, 3], [1, 3], [2, 3], [2, 2], [2, 1], [2, 0], [1, 0]])
