from sunpy.util import util


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


def test_expand_list_generator():
    lst = ['a', 'b', [], (['c', 'd']), tuple(), ['e']]
    assert list(util.expand_list_generator(lst)) == ['a', 'b', 'c', 'd', 'e']


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
