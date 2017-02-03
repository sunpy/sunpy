import pytest

from sunpy.util.metadata import MetaDict


class Base:

    def setup_method(self):

        self.seas = [['labrador', 'americas'],
                     ['Norwegian', 'Europe'],
                     ['BALTIC', 'europe'],
                     ['LaPteV', 'arctic']]

        # distinct keys to self.seas
        self.more_seas = [['Aegean', 'greece'],
                          ['CAtaLaN', 'mediterranean'],
                          ['IONIAN', 'italy'],
                          ['amundsen', 's.ocean']]

        # some similar, 'like',  keys to self.seas
        self.yet_more_seas = [['norwegian', 'NORWAY'],
                              ['SCOTIA', 's.ocean'],
                              ['laptev', 'russia'],
                              ['Mawson', 'E.Antartica']]

        self.md = MetaDict(self.seas)

    @staticmethod
    def check_contents(metadict_inst, expected):
        assert len(metadict_inst) == len(expected)

        for key, val in expected:
            assert metadict_inst[key.upper()] == val

    @staticmethod
    def check_insertion_order(metadict_inst, expected):

        def normalise_keys(lst_pairs):
            return [key.lower() for key, _ in lst_pairs]

        assert normalise_keys(metadict_inst.items()) == \
            normalise_keys(expected)

    @staticmethod
    def check(metadict_inst, expected):
        """ Check a MetaDict instance

        actual against expected contents and ensure insertion order is
        preserved.

        :param metadict_inst: a MetaDict
        :param expected: iterable of key/value pairs
        """
        Base.check_contents(metadict_inst, expected)
        Base.check_insertion_order(metadict_inst, expected)

    @staticmethod
    def pairs_to_dict(lst_of_pairs):
        """ convert a list/tuple of lists/tuples to a dictionary

        E.g.  [['a', 1], ['b', 2]] -> {'a': 1, 'b': 2}
        """
        return {kv_pair[0]: kv_pair[1] for kv_pair in lst_of_pairs}


class TestInitialisation(Base):
    """ All forms of instantiating a MetaDict
    """

    def test_with_dict(self):
        in_dict = TestInitialisation.pairs_to_dict(self.seas)

        # Using traditional dictionary, order of insertion is *not* preserved
        TestInitialisation.check_contents(MetaDict(in_dict), self.seas)

    def test_with_list_of_lists(self):
        TestInitialisation.check(MetaDict(self.more_seas), self.more_seas)

    def test_init_with_tuple_of_tuples(self):
        in_tuple = tuple(tuple(kv_pair) for kv_pair in self.seas)

        TestInitialisation.check(MetaDict(in_tuple), self.seas)

    def test_with_metadict(self):
        TestInitialisation.check(MetaDict(MetaDict(self.yet_more_seas)),
                                 self.yet_more_seas)

    def test_with_illegal_arg(self):
        with pytest.raises(TypeError):
            MetaDict(set(('a', 'b', 'c', 'd')))


class TestMethods(Base):

    def test_membership(self):
        assert 'LABRADOR' in self.md
        assert 'labrador' in self.md
        assert 'baLtiC' in self.md

        assert 'pechora' not in self.md

    def test_getitem(self):
        assert self.md['norwegian'] == 'Europe'
        assert self.md['NORWEGIAN'] != 'europe'
        assert self.md['laptev'] == 'arctic'

        assert self.md['Norwegian'] != 'europe'

        with pytest.raises(KeyError):
            self.md['key-not-exists']

    def test_get(self):
        assert self.md.get('baltic') == 'europe'
        assert self.md.get('atlantic') is None

        default_value = 'Unknown'
        assert self.md.get('sargasso', default=default_value) == default_value

    def test_setitem_existing(self):
        len_before = len(self.md)
        self.md['NORWEGIAN'] = 'Scandinavia'
        assert len(self.md) == len_before

        assert self.md['norwegian'] == 'Scandinavia'
        assert self.md['NoRwEgIaN'] == 'Scandinavia'

    def test_setitem_new(self):
        len_before = len(self.md)
        self.md['Irish'] = 'N.Europe'
        assert len(self.md) == len_before + 1

        assert self.md['Irish'] == 'N.Europe'
        assert self.md['irish'] == 'N.Europe'
        assert self.md['IRISH'] == 'N.Europe'

        assert self.md['IrisH'] != 'n.europe'

    def test_setdefault(self):
        self.md.setdefault('poseidon', 'NOT-FOUND')
        assert self.md.get('Poseidon') == 'NOT-FOUND'
        assert self.md['pOSeidon'] == 'NOT-FOUND'

        # Should only impact key 'poseidon'
        assert self.md.get('Wandel') is None

    def test_has_key(self):

        # suppress pep8/flake8 'W601 .has_key() is deprecated, ...' warnings
        # as MetaDict explicitly supports the 'has_key()' method
        assert self.md.has_key('laptev') is True  # noqa
        assert self.md.has_key('BALtIC') is True  # noqa

        assert self.md.has_key('Beaufort') is False  # noqa

    def test_pop(self):
        len_before = len(self.md)

        self.md.pop('kara') is None
        assert len(self.md) == len_before

        default_value = 'not-recognized'
        assert self.md.pop('lincoln', default=default_value) == default_value
        assert len(self.md) == len_before

        assert self.md.pop('baltic') == 'europe'
        assert len(self.md) == len_before - 1


class TestUpdateMethod(Base):
    """ MetaDict.update(...)
    """

    def setup(self):

        # Expected when self.seas & self.yet_more_seas are added to a MetaDict
        self.combined = [['labrador', 'americas'],
                         ['Norwegian', 'NORWAY'],
                         ['BALTIC', 'europe'],
                         ['LaPteV', 'russia'],
                         ['SCOTIA', 's.ocean'],
                         ['Mawson', 'E.Antartica']]

    def test_update_with_distinct_keys(self):
        self.md.update(MetaDict(self.more_seas))

        TestUpdateMethod.check(self.md, self.seas + self.more_seas)

    def test_update_with_like_keys(self):
        # values of existing keys should be updated but, original insertion
        # order should be preserved
        self.md.update(MetaDict(self.yet_more_seas))

        TestUpdateMethod.check(self.md, self.combined)

    def test_update_with_dict(self):
        # Updating with a dictionary, order is *not* preserved
        self.md.update(TestUpdateMethod.pairs_to_dict(self.yet_more_seas))

        TestUpdateMethod.check_contents(self.md, self.combined)

    def test_update_with_same_metadict(self):
        self.md.update(self.md)

        TestUpdateMethod.check(self.md, self.seas)
