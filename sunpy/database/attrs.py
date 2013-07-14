from sunpy.net.attr import AttrWalker, Attr, ValueAttr, AttrAnd, AttrOr
from sunpy.database import tables

__all__ = ['Starred', 'Tag', 'walker']


class _BooleanAttr(object):
    def __init__(self, value, make):
        self.value = bool(value)
        self.make = make

    def __and__(self, other):
        if not isinstance(other, self.make):
            return AttrAnd([self, other])
        attr = self.make()
        attr.value = self.value and other.value
        return attr

    def __or__(self, other):
        if not isinstance(other, self.make):
            return AttrOr([self, other])
        attr = self.make()
        attr.value = self.value or other.value
        return attr

    def __nonzero__(self):
        return self.value

    def __invert__(self):
        attr = self.make()
        attr.value = not self.value
        return attr

    def __eq__(self, other):
        return isinstance(other, self.make) and self.value == other.value

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<%s(%s)>' % (self.__class__.__name__, self.value)


class Starred(_BooleanAttr, Attr):
    def __init__(self):
        _BooleanAttr.__init__(self, True, self.__class__)


class Tag(Attr):
    def __init__(self, tagname):
        self.tagname = tagname
        self.inverted = False

    def __invert__(self):
        tag = self.__class__(self.tagname)
        tag.inverted = True
        return tag

    def collides(self, other):  # pragma: no cover
        return False

    def __repr__(self):
        return '<Tag(%r, %r)>' % (self.tagname, self.inverted)


walker = AttrWalker()


@walker.add_creator(AttrOr)
def _create(wlk, root, session):
    entries = set()
    for attr in root.attrs:
        # make sure to add only new entries to the set to avoid duplicates
        entries.update(set(wlk.create(attr, session)) - entries)
    return list(entries)


@walker.add_creator(AttrAnd)
def _create(wlk, root, session):
    entries = [set(wlk.create(attr, session)) for attr in root.attrs]
    return list(set.intersection(*entries))


@walker.add_creator(ValueAttr)
def _create(wlk, root, session):
    query = session.query(tables.DatabaseEntry)
    for key, value in root.attrs.iteritems():
        typ = key[0]
        if typ == 'tag':
            criterion = tables.Tag.name.in_([value])
            # `key[1]` is here the `inverted` attribute of the tag. That means
            # that if `subtyp` is True, the given tag must not be included in
            # the resulting entries.
            if key[1]:
                query = query.filter(~tables.DatabaseEntry.tags.any(criterion))
            else:
                query = query.filter(tables.DatabaseEntry.tags.any(criterion))
        else:
            query = query.filter_by(**{typ: value})
    return query.all()


@walker.add_converter(Tag)
def _convert(attr):
    return ValueAttr({('tag', attr.inverted): attr.tagname})


@walker.add_converter(Starred)
def _convert(attr):
    return ValueAttr({('starred', ): attr.value})
