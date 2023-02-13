
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr

walker = AttrWalker()


@classmethod
def register_values(cls):

    from sunpy.net import attrs
    adict = {
        attrs.Instrument: [('DSLP', 'Dual Segmented Langmuir Probe')],
        attrs.Source: [('ESA', 'European Space Agency')],
        # attrs.Provider: [('SDAC', 'Solar Data Analysis Center')],
        # attrs.Detector: [('C1', 'Coronograph 1'),
        #                  ('C2', 'Coronograph 2'),
        #                  ('C3', 'Coronograph 3')]
    }

    return adict


@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    results = []
    for sub in tree.attrs:
        results.append(wlk.create(sub))

    return results


@walker.add_creator(AttrAnd, DataAttr)
def create_and(wlk, tree):
    result = dict()
    wlk.apply(tree, result)
    return [result]

# @walker.add_applier(a.Time)
# def _(wlk, attr, params):
#     return params.update({'startTime': attr.start.isot,
#                           'endTime': attr.end.isot})

# @walker.add_applier(a.Level)
# def _(wlk, attr, params):
#     return params.update({'level': attr.value})
