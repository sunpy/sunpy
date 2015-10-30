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
