import copy

import pytest

from sunpy.util.metadata import MetaDict, ModifiedItem


# TODO: Should these be here and not in util?
def check_contents(metadict_inst, expected):
    """
    Ensure that the key/values of `metadict_inst` match those of the key/value
    pairs in `expected`.

    The case of the keys is ignored, as is the order.
    """
    assert len(metadict_inst) == len(expected)

    for key, val in expected:
        assert metadict_inst[key.upper()] == val


def check_insertion_order(metadict_inst, expected):
    """
    Ensure that the keys of `metadict_inst` are in the same order as the keys
    in `expected`.

    The keys case is ignored.
    """
    def normalise_keys(lst_pairs):
        return [key.lower() for key, _ in lst_pairs]

    assert normalise_keys(metadict_inst.items()) == normalise_keys(expected)


def check_contents_and_insertion_order(metadict_inst, expected):
    """
    Ensure both the contents and order of the the key/values of `metadict_inst`
    match the key/values in `expected`. Keys are case insensitive.

    Parameters
    ----------
    metadict_inst : sunpy.util.metadata.MetaDict
        a Metadict instance under test
    expected : iterable object of key/value pairs
        the values we expect from a test
    """
    check_contents(metadict_inst, expected)
    check_insertion_order(metadict_inst, expected)


def pairs_to_dict(lst_of_pairs):
    """
    Convert a list/tuple of lists/tuples to a dictionary.

    Parameters
    ----------
    lst_of_pairs : iteterable object of iterable objects
        an iterable containing iterables, each of these
        contained iterables is a key/value pair.

    Examples
    --------
    >>> pairs_to_dict([['a', 1], ['b', 2]])   # doctest: +SKIP
    {'a': 1, 'b': 2}
    >>> pairs_to_dict([('x', 34), ('y', 56)])   # doctest: +SKIP
    {'x': 34, 'y': 56}
    """
    return {kv_pair[0]: kv_pair[1] for kv_pair in lst_of_pairs}


@pytest.fixture
def sea_locations():
    return [['labrador', 'americas'],
            ['Norwegian', 'Europe'],
            ['BALTIC', 'europe'],
            ['LaPteV', 'arctic']]


@pytest.fixture
def seas_metadict(sea_locations):
    return MetaDict(sea_locations)


@pytest.fixture
def atomic_weights():
    return [['hydrogen', 1],
            ['chromium', 24],
            ['mercury', 80],
            ['iridium', 77]]


@pytest.fixture
def atomic_weights_keycomments():
    """
    Atomic weights with keycomments (including extraneous keycomments).
    """
    return [['hydrogen', 1],
            ['chromium', 24],
            ['mercury', 80],
            ['iridium', 77],
            ['keycomments', {
                'chromium': 'Cr',
                'extra key 1': 'foo',
                'MERCURY': 'Hg',
                'extra key 2': 'bar'
            }]]


@pytest.fixture
def atomic_weights_pruned_keycomments():
    """
    Expected result of atomic weights with extraneous keycomments removed.
    """
    return [['hydrogen', 1],
            ['chromium', 24],
            ['mercury', 80],
            ['iridium', 77],
            ['keycomments', {
                'chromium': 'Cr',
                'MERCURY': 'Hg'
            }]]


@pytest.fixture
def atomic_weights_keycomments_after_removal():
    """
    Expected result after removing some atomic weights and their keycomments.
    """
    return [['chromium', 24],
            ['iridium', 77],
            ['keycomments', {
                'chromium': 'Cr'
            }]]


@pytest.fixture
def empty_keycomments():
    """
    Expected result after removing all keys.
    """
    return [['keycomments', {}]]


# Test constructors `MetaDict.__init__(...)`


def test_init_with_lst_of_lsts(seas_metadict, sea_locations):
    """
    Initialise "Metadict" with a list of lists.

    Each sub list is a key/value pair.

    Examples
    --------
    >>> m = MetaDict([['H', 'Hydrogen'], ['B', 'Boron']])
    """
    check_contents_and_insertion_order(MetaDict(sea_locations), sea_locations)


def test_init_with_tuple_of_tuples(sea_locations):
    """
    Initialise "Metadict" with a tuple of tuples.

    Each 'sub-tuple' is a key/value pair.

    Examples
    --------
    >>> c = MetaDict((('Be', 21), ('Nb', 21)))
    """
    in_tuple = tuple(tuple(kv_pair) for kv_pair in sea_locations)
    check_contents_and_insertion_order(MetaDict(in_tuple), sea_locations)


def test_init_with_dict(sea_locations):
    """
    Initialise "Metadict" with standard Python dictionary.

    Order of insertion is *not* preserved - only check the contents.
    """
    check_contents(MetaDict(pairs_to_dict(sea_locations)), sea_locations)


def test_init_with_metadict(atomic_weights):
    """
    Initialise "Metadict" with another "Metadict".
    """
    original = MetaDict(atomic_weights)
    new = MetaDict(original)
    check_contents_and_insertion_order(new, atomic_weights)

    # It's not just a simple reassignemnet
    assert id(original) != id(new)


def test_init_with_keycomments(atomic_weights_keycomments, atomic_weights_pruned_keycomments):
    """
    Initialise "Metadict" with keycomments. Ensure caller's keycomments dict is not mutated.
    """
    orig_dict = pairs_to_dict(atomic_weights_keycomments)
    orig_keycomments = orig_dict['keycomments'].copy()

    md = MetaDict(orig_dict)
    check_contents_and_insertion_order(md, atomic_weights_pruned_keycomments)

    assert md['keycomments'] is not orig_dict['keycomments']
    assert orig_dict['keycomments'] == orig_keycomments


def test_init_with_illegal_arg():
    """
    Ensure attempt to initialise with a nonsensical data structure is rejected.
    """
    with pytest.raises(TypeError):
        MetaDict({'a', 'b', 'c', 'd'})


def test_init_with_invalid_keycomments_type():
    """
    Ensure attempt to initialise with an invalid keycomments type is rejected.
    """
    with pytest.raises(TypeError):
        MetaDict({'a': 1, 'b': 2, 'keycomments': 3})


# Test individual methods


def test_membership(seas_metadict):
    """
    Test `MetaDict.__contains__(...)`
    """
    assert 'labrador' in seas_metadict
    assert 'pechora' not in seas_metadict


def test_getitem(seas_metadict):
    """
    Test `MetaDict.__getitem__(...)`
    """
    assert seas_metadict['Norwegian'] == 'Europe'
    assert seas_metadict['BALTIC'] != 'arctic'

    with pytest.raises(KeyError):
        seas_metadict['key-not-exists']


def test_setitem_existing_entry(seas_metadict):
    """
    Test `MetaDict.__setitem__(...)` on an existing entry.
    """
    len_before = len(seas_metadict)
    seas_metadict['NORWEGIAN'] = 'Scandinavia'
    assert len(seas_metadict) == len_before

    assert seas_metadict['NORWEGIAN'] == 'Scandinavia'


def test_setitem_new_entry(seas_metadict):
    """
    Test `MetaDict.__setitem__(...)`.

    Add a new entry
    """
    len_before = len(seas_metadict)
    seas_metadict['Irish'] = 'N.Europe'
    assert len(seas_metadict) == len_before + 1

    assert seas_metadict['Irish'] == 'N.Europe'


def test_delitem(seas_metadict):
    """
    Test `MetaDict.__delitem__(...)`.
    """
    len_before = len(seas_metadict)
    del seas_metadict['NoRwEgIaN']
    del seas_metadict['baltic']
    assert len(seas_metadict) == len_before - 2

    with pytest.raises(KeyError):
        seas_metadict['baltic']

    with pytest.raises(KeyError):
        seas_metadict['NoRwEgIaN']


def test_delitem_missing_key(seas_metadict):
    """
    Test `MetaDict.__delitem__(...)` raises error on missing key.
    """
    with pytest.raises(KeyError):
        del seas_metadict['missing key']


def test_delitem_with_keycomments(atomic_weights_keycomments,
                                  atomic_weights_keycomments_after_removal):
    """
    Test `MetaDict.__delitem__(...)` removes corresponding keycomments.
    """
    len_before = len(atomic_weights_keycomments)
    md = MetaDict(atomic_weights_keycomments)
    del md['hydrogen']
    del md['mercury']
    assert len(md) == len_before - 2

    with pytest.raises(KeyError):
        md['hydrogen']

    with pytest.raises(KeyError):
        md['mercury']

    check_contents_and_insertion_order(md, atomic_weights_keycomments_after_removal)


def test_get(seas_metadict):
    """
    Test `MetaDict.get(...)`
    """
    assert seas_metadict.get('BALTIC') == 'europe'
    assert seas_metadict.get('atlantic') is None

    default_value = 'Unknown'
    assert seas_metadict.get('sargasso', default=default_value) == default_value


def test_setdefault(seas_metadict):
    """
    Test `MetaDict.setdefault(...)`
    """
    seas_metadict.setdefault('poseidon', 'NOT-FOUND')
    assert seas_metadict.get('poseidon') == 'NOT-FOUND'
    assert seas_metadict['poseidon'] == 'NOT-FOUND'

    # Should only impact key 'poseidon'
    assert seas_metadict.get('Wandel') is None


def test_has_key(seas_metadict):
    """
    Test `MetaDict.has_key(...)`
    """
    # MetaDict explicitly supports the 'has_key()' method
    assert seas_metadict.has_key('LaPteV') is True
    assert seas_metadict.has_key('Beaufort') is False


def test_pop(seas_metadict):
    """
    Test `MetaDict.pop(...)`
    """
    # Nothing to 'pop', nothing should change
    len_before = len(seas_metadict)
    seas_metadict.pop('kara') is None
    assert len(seas_metadict) == len_before

    # Nothing to 'pop', nothing should change but, we should get the default value
    default_value = 'not-recognized'
    assert seas_metadict.pop('lincoln', default=default_value) == default_value
    assert len(seas_metadict) == len_before

    assert seas_metadict.pop('baltic') == 'europe'
    assert len(seas_metadict) == len_before - 1


def test_pop_with_keycomments(atomic_weights_keycomments,
                              atomic_weights_keycomments_after_removal):
    """
    Test `MetaDict.pop(...)` removes corresponding keycomments.
    """
    len_before = len(atomic_weights_keycomments)
    md = MetaDict(atomic_weights_keycomments)

    assert md.pop('hydrogen') == 1
    assert md.pop('mercury') == 80
    assert len(md) == len_before - 2
    check_contents_and_insertion_order(md, atomic_weights_keycomments_after_removal)


def test_popitem_with_keycomments(atomic_weights_keycomments, empty_keycomments):
    """
    Test `MetaDict.popitem(...)` removes corresponding keycomments.
    """
    md = MetaDict(atomic_weights_keycomments)

    assert md.popitem(last=False) == ('hydrogen', 1)
    assert md.popitem(last=False) == ('chromium', 24)
    assert md.popitem(last=False) == ('mercury', 80)
    assert md.popitem(last=False) == ('iridium', 77)
    check_contents_and_insertion_order(md, empty_keycomments)

#  Test `MetaDict.update(...)`.


def test_update_with_distinct_keys(seas_metadict, sea_locations, atomic_weights):
    """
    Update the MetaDict 'world_seas', with another MetaDict, 'chem_elems' which
    has a completely different set of keys.

    This should result in adding 'chem_elems' onto 'world_seas'.
    Original insertion order of 'world_seas' should be preserved.
    """
    world_seas = seas_metadict
    chem_elems = MetaDict(atomic_weights)
    world_seas.update(chem_elems)

    check_contents_and_insertion_order(world_seas, sea_locations + atomic_weights)


@pytest.fixture
def seas_and_atomic_weights():
    """
    A mixture of key/value data from `seas` & `atomic_weights`
    """
    return [['norwegian', 'NORWAY'],
            ['Chromium', 24],
            ['laptev', 'russia'],
            ['Iridium', 77]]


@pytest.fixture
def combined_seas_atomic():
    """
    The expected result of a "Metadict" initialized with `sea_locations` and
    then updated with `seas_and_atomic_weights`
    """
    return [['labrador', 'americas'],
            ['Norwegian', 'NORWAY'],
            ['BALTIC', 'europe'],
            ['LaPteV', 'russia'],
            ['Chromium', 24],
            ['Iridium', 77]]


def test_update_with_like_keys(seas_metadict, seas_and_atomic_weights, combined_seas_atomic):
    """
    Update the "Metadict" 'world_seas' with another "Metadict", 'atomic_seas'.

    'atomic_seas' has some keys which are the same as 'world_seas', some
    are different.

    In 'world_seas', values of existing keys should be updated but, original
    insertion order should be preserved
    """
    world_seas = seas_metadict
    atomic_seas = MetaDict(seas_and_atomic_weights)
    world_seas.update(atomic_seas)
    check_contents_and_insertion_order(world_seas, combined_seas_atomic)


def test_update_with_dict(seas_metadict, seas_and_atomic_weights, combined_seas_atomic):
    """
    Update an existing "Metadict" with a standard python dictionary some of
    whose keys are the same, some are different.

    In the updated "Metadict", values of existing keys should be updated
    but as we are using a standard dictionary, insertion order of the
    new items is non-deterministic so only check the contents of the
    updated structure.
    """
    seas_metadict.update(pairs_to_dict(seas_and_atomic_weights))
    check_contents(seas_metadict, combined_seas_atomic)


def test_update_with_same_metadict(seas_metadict, sea_locations):
    """
    Update a 'MetaDict' with itself.

    Nothing should changes.
    """
    seas_metadict.update(seas_metadict)

    check_contents_and_insertion_order(seas_metadict, sea_locations)


def test_key_case_insensitivity(seas_metadict):
    """
    The keys of a "Metadict" are case insensitive.

    Using the key 'BALTIC' is identical to the key 'baltic'.
    """
    # membership
    assert 'laptev' in seas_metadict
    assert 'LAPTEV' in seas_metadict

    # get
    assert seas_metadict['norwegian'] == 'Europe'
    assert seas_metadict['NORWEGIAN'] == 'Europe'

    assert seas_metadict.get('labrador') == 'americas'
    assert seas_metadict.get('labRAdor') == 'americas'

    # MetaDict explicitly supports the 'has_key()' method
    assert seas_metadict.has_key('BALTIC')
    assert seas_metadict.has_key('balTIC')

    # This key already exists. It should *not* take on the default value
    seas_metadict.setdefault('norwEgiaN', default='Sweden')
    assert seas_metadict['norwEgiaN'] == 'Europe'
    assert seas_metadict['Norwegian'] == 'Europe'
    assert seas_metadict.get('norwEgiaN') == 'Europe'

    # setting, existing entry
    seas_metadict['LaPteV'] = 'japan'
    assert seas_metadict['LaPteV'] == 'japan'
    assert seas_metadict['LAPTEV'] == 'japan'
    assert seas_metadict.get('LApTEV') == 'japan'

    # setting, new entry
    seas_metadict['bering'] = 'Russia'
    assert seas_metadict['bering'] == 'Russia'
    assert seas_metadict['BeRinG'] == 'Russia'
    assert seas_metadict.get('BERING') == 'Russia'


def test_original_copy():
    md = MetaDict({'foo': 'bar'})
    assert md.original_meta == md

    # Add a key, make sure original contents doesn't change
    md['a'] = 'b'
    assert md.original_meta != md
    assert list(md.keys()) == ['foo', 'a']
    assert list(md.original_meta.keys()) == ['foo']
    assert md.added_items == {'a': 'b'}

    # Check that creating a new MetaDict preserves the original copy
    md = MetaDict(md)
    assert list(md.keys()) == ['foo', 'a']
    assert list(md.original_meta.keys()) == ['foo']

    # Check that creating a copy preserves the original copy
    md = copy.copy(md)
    assert list(md.keys()) == ['foo', 'a']
    assert list(md.original_meta.keys()) == ['foo']

    # Check that creating a deepcopy preserves the original copy
    md = copy.deepcopy(md)
    assert list(md.keys()) == ['foo', 'a']
    assert list(md.original_meta.keys()) == ['foo']

    # Check that creating using .copy() preserves the original copy
    md = md.copy()
    assert list(md.keys()) == ['foo', 'a']
    assert list(md.original_meta.keys()) == ['foo']

    # Check modification of items
    md['foo'] = 'bar1'
    assert md.modified_items == {'foo': ModifiedItem('bar', 'bar1')}

    # Check removal of items
    md.pop('foo')
    assert md.removed_items == {'foo': 'bar'}
