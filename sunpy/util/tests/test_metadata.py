import pytest

from sunpy.util.metadata import MetaDict


def check_contents(metadict_inst, expected):
    assert len(metadict_inst) == len(expected)

    for key, val in expected:
        assert metadict_inst[key.upper()] == val


def check_insertion_order(metadict_inst, expected):

    def normalise_keys(lst_pairs):
        return [key.lower() for key, _ in lst_pairs]

    assert normalise_keys(metadict_inst.items()) == normalise_keys(expected)


def check(metadict_inst, expected):
    """
    Check a MetaDict instance

    Ensure content and insertion order are preserved

    Parameters
    ----------
    metadict_inst: sunpy.util.metadata.MetaDict
                  a Metadict instance under test

    expected: iterable object of key/value pairs
             the value we expect from a test
    """
    check_contents(metadict_inst, expected)
    check_insertion_order(metadict_inst, expected)


def pairs_to_dict(lst_of_pairs):
    """
    Convert a list/tuple of lists/tuples to a dictionary

    Parameters
    ----------
    lst_of_pairs: iteterable object of iterable objects
                  an iterable containing iterables, each of these
                  contained iterables is a key/value pair.

    Examples
    --------
    >>> pairs_to_dict([['a', 1], ['b', 2]])
    {'a': 1, 'b': 2}
    >>> pairs_to_dict([('x', 34), ('y', 56)])
    {'x': 34, 'y': 56}
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
    Contains keys distinct from 'std_test_data'
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


def test_init_with_lst_of_lsts(std_test_data):
    check(MetaDict(std_test_data), std_test_data)


def test_init_with_tuple_of_tuples(std_test_data):
    in_tuple = tuple(tuple(kv_pair) for kv_pair in std_test_data)
    check(MetaDict(in_tuple), std_test_data)


def test_init_with_dict(std_test_data):
    # Using traditional dictionary, order of insertion is *not* preserved
    check_contents(MetaDict(pairs_to_dict(std_test_data)), std_test_data)


def test_init_with_metadict(same_diff_keys_data):
    check(MetaDict(MetaDict(same_diff_keys_data)), same_diff_keys_data)


def test_init_with_illegal_arg():
    with pytest.raises(TypeError):
        MetaDict(set(('a', 'b', 'c', 'd')))


def test_membership_op(std_metadict):
    assert 'LABRADOR' in std_metadict
    assert 'labrador' in std_metadict
    assert 'baLtiC' in std_metadict

    assert 'pechora' not in std_metadict


def test_getitem_op(std_metadict):
    assert std_metadict['norwegian'] == 'Europe'
    assert std_metadict['NORWEGIAN'] != 'europe'
    assert std_metadict['laptev'] == 'arctic'

    assert std_metadict['Norwegian'] != 'europe'

    with pytest.raises(KeyError):
        std_metadict['key-not-exists']


def test_get_method(std_metadict):
    assert std_metadict.get('baltic') == 'europe'
    assert std_metadict.get('atlantic') is None

    default_value = 'Unknown'
    assert std_metadict.get('sargasso', default=default_value) == default_value


def test_setitem_op_existing(std_metadict):
    """
    Update an existing entry
    """
    len_before = len(std_metadict)
    std_metadict['NORWEGIAN'] = 'Scandinavia'
    assert len(std_metadict) == len_before

    assert std_metadict['norwegian'] == 'Scandinavia'
    assert std_metadict['NoRwEgIaN'] == 'Scandinavia'


def test_setitem_op_new(std_metadict):
    """
    Add a new entry
    """
    len_before = len(std_metadict)
    std_metadict['Irish'] = 'N.Europe'
    assert len(std_metadict) == len_before + 1

    assert std_metadict['Irish'] == 'N.Europe'
    assert std_metadict['irish'] == 'N.Europe'
    assert std_metadict['IRISH'] == 'N.Europe'

    assert std_metadict['IrisH'] != 'n.europe'


def test_setdefault_method(std_metadict):
    std_metadict.setdefault('poseidon', 'NOT-FOUND')
    assert std_metadict.get('Poseidon') == 'NOT-FOUND'
    assert std_metadict['pOSeidon'] == 'NOT-FOUND'

    # Should only impact key 'poseidon'
    assert std_metadict.get('Wandel') is None


def test_has_key_method(std_metadict):

    # suppress pep8/flake8 'W601 .has_key() is deprecated, ...' warnings
    # as MetaDict explicitly supports the 'has_key()' method
    assert std_metadict.has_key('laptev') is True  # noqa
    assert std_metadict.has_key('BALtIC') is True  # noqa

    assert std_metadict.has_key('Beaufort') is False  # noqa


def test_pop_method(std_metadict):
    len_before = len(std_metadict)

    std_metadict.pop('kara') is None
    assert len(std_metadict) == len_before

    default_value = 'not-recognized'
    assert std_metadict.pop('lincoln', default=default_value) == default_value
    assert len(std_metadict) == len_before

    assert std_metadict.pop('baltic') == 'europe'
    assert len(std_metadict) == len_before - 1


def test_update_method_with_distinct_keys(std_metadict, std_test_data,
                                          differet_keys_data):

    std_metadict.update(MetaDict(differet_keys_data))

    check(std_metadict, std_test_data + differet_keys_data)


def test_update_method_with_like_keys(std_metadict, same_diff_keys_data,
                                      combined_data):
    # values of existing keys should be updated but, original insertion
    # order should be preserved

    std_metadict.update(MetaDict(same_diff_keys_data))

    check(std_metadict, combined_data)


def test_update_method_with_dict(std_metadict, same_diff_keys_data,
                                 combined_data):
    # Updating with a dictionary, order is *not* preserved
    std_metadict.update(pairs_to_dict(same_diff_keys_data))

    check_contents(std_metadict, combined_data)


def test_update_method_with_same_metadict(std_metadict, std_test_data):
    std_metadict.update(std_metadict)

    check(std_metadict, std_test_data)
