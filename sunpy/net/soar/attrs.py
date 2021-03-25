import sunpy.net.attrs as a
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, SimpleAttr


class Identifier(SimpleAttr):
    """
    The data product identifier to search for.
    """


walker = AttrWalker()


@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    results = []
    for sub in tree.attrs:
        results.append(wlk.create(sub))
    return results


@walker.add_creator(AttrAnd)
def create_and(wlk, tree):
    result = []
    wlk.apply(tree, result)
    return [result]


@walker.add_applier(AttrAnd)
def apply_and(wlk, and_attr, params):
    for iattr in and_attr.attrs:
        wlk.apply(iattr, params)


@walker.add_applier(a.Time)
def _(wlk, attr, params):
    start = attr.start.strftime('%Y-%m-%d+%H:%M:%S')
    end = attr.end.strftime('%Y-%m-%d+%H:%M:%S')
    params.append(f"begin_time>='{start}'+AND+begin_time<='{end}'")


@walker.add_applier(a.Level)
def _(wlk, attr, params):
    level = int(attr.value)
    valid_levels = [0, 1, 2, 3]
    if int(level) not in valid_levels:
        raise ValueError(f'Valid Solar Orbiter data levels are {valid_levels}')
    params.append(f"level='L{level}'")


@walker.add_applier(a.Instrument)
def _(wlk, attr, params):
    params.append(f"instrument='{attr.value}'")


@walker.add_applier(Identifier)
def _(wlk, attr, params):
    params.append(f"descriptor='{attr.value}'")
