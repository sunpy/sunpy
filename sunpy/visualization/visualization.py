# -*- coding: utf-8 -*-
from matplotlib import pyplot


def toggle_pylab(fn):
    """
    A decorator to prevent functions from opening matplotlib windows
    unexpectedly when sunpy is run in interactive shells like ipython
    --pylab.

    Toggles the value of matplotlib.pyplot.isinteractive() to preserve the
    users' expectations of pylab's behaviour in general.

    Parameters
    ----------
    fn : function object
        ?

    Returns
    ------
    ? : ?
        ?
    .. todo::
        improve documentation
    """

    if pyplot.isinteractive():
        def fn_itoggle(*args, **kwargs):
            pyplot.ioff()
            ret = fn(*args, **kwargs)
            pyplot.ion()
            return ret
        return fn_itoggle
    else:
        return fn


def axis_labels_from_ctype(ctype, unit):
    ctype_short = ctype[:4]

    labels = {'HGLN': 'Heliographic Longitude [{}]'.format(unit),
              'CRLN': 'Carrington Longitude [{}]'.format(unit),
              'HPLN': 'Helioprojective Longitude (Solar-X) [{}]'.format(unit),
              'SOLX': 'Heliocentric X [{}]'.format(unit),

              'HGLT': 'Latitude [{}]'.format(unit),
              'CRLT': 'Latitude [{}]'.format(unit),
              'HPLT': 'Helioprojective Latitude (Solar-Y) [{}]'.format(unit),
              'SOLY': 'Heliocentric Y [{}]'.format(unit)}

    return labels.get(ctype_short, "{} [{}]".format(ctype, unit))
