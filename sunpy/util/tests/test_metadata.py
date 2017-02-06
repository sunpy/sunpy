import pytest

from sunpy.util.metadata import MetaDict


def check_contents(metadict_inst, expected):
    assert len(metadict_inst) == len(expected)

    for key, val in expected:
        assert metadict_inst[key.upper()] == val


def check_insertion_order(metadict_inst, expected):

    def normalise_keys(lst_pairs):
        return [key.lower() for key, _ in lst_pairs]

    assert normalise_keys(metadict_inst.items()) == \
        normalise_keys(expected)


def check(metadict_inst, expected):
    """ Check a MetaDict instance

    actual against expected contents and ensure insertion order is
    preserved.

    :param metadict_inst: a MetaDict
    :param expected: iterable of key/value pairs
    """
    check_contents(metadict_inst, expected)
    check_insertion_order(metadict_inst, expected)


def pairs_to_dict(lst_of_pairs):
    """ convert a list/tuple of lists/tuples to a dictionary

    E.g.  [['a', 1], ['b', 2]] -> {'a': 1, 'b': 2}
    """
    return {kv_pair[0]: kv_pair[1] for kv_pair in lst_of_pairs}


@pytest.fixture
def std_test_data():
    return [['labrador', 'americas'],
            ['Norwegian', 'Europe'],
            ['BALTIC', 'europe'],
            ['LaPteV', 'arctic']]


@pytest.fixture
def differet_keys_data():
    """
    Contains distinct keys from 'std_test_data'
    """
    return [['Aegean', 'greece'],
            ['CAtaLaN', 'mediterranean'],
            ['IONIAN', 'italy'],
            ['amundsen', 's.ocean']]


@pytest.fixture
def same_diff_keys_data():
    """
    Contains a mixture of same and distinct keys from 'std_test_data'
    """
    return [['norwegian', 'NORWAY'],
            ['SCOTIA', 's.ocean'],
            ['laptev', 'russia'],
            ['Mawson', 'E.Antartica']]


@pytest.fixture
def combined_data():
    """
    The expected result of combining the 'std_test_data' &
    'same_diff_keys_data' in a MetaDict
    """
    return [['labrador', 'americas'],
            ['Norwegian', 'NORWAY'],
            ['BALTIC', 'europe'],
            ['LaPteV', 'russia'],
            ['SCOTIA', 's.ocean'],
            ['Mawson', 'E.Antartica']]


@pytest.fixture
def std_metadict(std_test_data):
    return MetaDict(std_test_data)


class TestInit:
    """
    All the ways on instantiating a MetaDict
    """

    def test_lst_of_lsts(self, std_test_data):
        check(MetaDict(std_test_data), std_test_data)

    def test_init_with_tuple_of_tuples(self, std_test_data):
        in_tuple = tuple(tuple(kv_pair) for kv_pair in std_test_data)
        check(MetaDict(in_tuple), std_test_data)

    def test_with_dict(self, std_test_data):
        # Using traditional dictionary, order of insertion is *not* preserved
        check_contents(MetaDict(pairs_to_dict(std_test_data)), std_test_data)

    def test_with_metadict(self, same_diff_keys_data):
        check(MetaDict(MetaDict(same_diff_keys_data)), same_diff_keys_data)

    def test_with_illegal_arg(self):
        with pytest.raises(TypeError):
            MetaDict(set(('a', 'b', 'c', 'd')))


class TestMethods:
    """
    Most of the API of MetaDict barring the initialization and 'update' methods
    """

    def test_membership(self, std_metadict):
        assert 'LABRADOR' in std_metadict
        assert 'labrador' in std_metadict
        assert 'baLtiC' in std_metadict

        assert 'pechora' not in std_metadict

    def test_getitem(self, std_metadict):
        assert std_metadict['norwegian'] == 'Europe'
        assert std_metadict['NORWEGIAN'] != 'europe'
        assert std_metadict['laptev'] == 'arctic'

        assert std_metadict['Norwegian'] != 'europe'

        with pytest.raises(KeyError):
            std_metadict['key-not-exists']

    def test_get(self, std_metadict):
        assert std_metadict.get('baltic') == 'europe'
        assert std_metadict.get('atlantic') is None

        default_value = 'Unknown'
        assert std_metadict.get('sargasso', default=default_value) \
            == default_value

    def test_setitem_existing(self, std_metadict):
        len_before = len(std_metadict)
        std_metadict['NORWEGIAN'] = 'Scandinavia'
        assert len(std_metadict) == len_before

        assert std_metadict['norwegian'] == 'Scandinavia'
        assert std_metadict['NoRwEgIaN'] == 'Scandinavia'

    def test_setitem_new(self, std_metadict):
        len_before = len(std_metadict)
        std_metadict['Irish'] = 'N.Europe'
        assert len(std_metadict) == len_before + 1

        assert std_metadict['Irish'] == 'N.Europe'
        assert std_metadict['irish'] == 'N.Europe'
        assert std_metadict['IRISH'] == 'N.Europe'

        assert std_metadict['IrisH'] != 'n.europe'

    def test_setdefault(self, std_metadict):
        std_metadict.setdefault('poseidon', 'NOT-FOUND')
        assert std_metadict.get('Poseidon') == 'NOT-FOUND'
        assert std_metadict['pOSeidon'] == 'NOT-FOUND'

        # Should only impact key 'poseidon'
        assert std_metadict.get('Wandel') is None

    def test_has_key(self, std_metadict):

        # suppress pep8/flake8 'W601 .has_key() is deprecated, ...' warnings
        # as MetaDict explicitly supports the 'has_key()' method
        assert std_metadict.has_key('laptev') is True  # noqa
        assert std_metadict.has_key('BALtIC') is True  # noqa

        assert std_metadict.has_key('Beaufort') is False  # noqa

    def test_pop(self, std_metadict):
        len_before = len(std_metadict)

        std_metadict.pop('kara') is None
        assert len(std_metadict) == len_before

        default_value = 'not-recognized'
        assert std_metadict.pop('lincoln', default=default_value) \
            == default_value
        assert len(std_metadict) == len_before

        assert std_metadict.pop('baltic') == 'europe'
        assert len(std_metadict) == len_before - 1


class TestUpdateMethod:

    def test_update_with_distinct_keys(self, std_metadict, std_test_data,
                                       differet_keys_data):

        std_metadict.update(MetaDict(differet_keys_data))

        check(std_metadict, std_test_data + differet_keys_data)

    def test_update_with_like_keys(self, std_metadict, same_diff_keys_data,
                                   combined_data):
        # values of existing keys should be updated but, original insertion
        # order should be preserved

        std_metadict.update(MetaDict(same_diff_keys_data))

        check(std_metadict, combined_data)

    def test_update_with_dict(self, std_metadict, same_diff_keys_data,
                              combined_data):
        # Updating with a dictionary, order is *not* preserved
        std_metadict.update(pairs_to_dict(same_diff_keys_data))

        check_contents(std_metadict, combined_data)

    def test_update_with_same_metadict(self, std_metadict, std_test_data):
        std_metadict.update(std_metadict)

        check(std_metadict, std_test_data)
