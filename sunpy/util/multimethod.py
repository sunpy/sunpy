"""
This module provides multimethod implementation in pure Python.
"""
from warnings import warn

from sunpy.util.exceptions import SunpyDeprecationWarning

__all__ = ['TypeWarning', 'MultiMethod']

SILENT = 0
WARN = 1
FAIL = 2


def _fmt_t(types):
    return ', '.join(type_.__name__ for type_ in types)


class TypeWarning(UserWarning):
    pass


class MultiMethod:
    """
    A multimethod is a callable object that decides which code to execute based
    on the type of one or more of its arguments.

    Parameters
    ----------
    get : `function`
        The function which receives args and kwargs and returns a tuple of values to consider for dispatch.
    """
    def __init__(self, get):
        warn("MultiMethod is deprecated and will be removed in sunpy 2.1", SunpyDeprecationWarning)
        self.get = get

        self.methods = []
        self.cache = {}

    def add(self, fun, types, override=SILENT):
        """
        Add ``fun`` to the multimethod. It will be executed if get returns
        values of the types passed as types. Must return tuples of same length
        for any input.

        Parameters
        ----------
        fun : `function`
            Function to be added to the multimethod.
        types : `tuple` of classes
            Types for which the function is executed.
        override : {SILENT, WARN, FAIL}
            Control behavior when overriding existing definitions.
            If it is set to ``SILENT``, prior definitions are silently
            overridden, if it is set to ``WARN`` a `sunpy.util.multimethod.TypeWarning`
            will be issued, and with ``FAIL`` a `TypeError` is raised when
            attempting to override an existing definition.
        """
        if override not in (SILENT, WARN, FAIL):
            raise ValueError(f"Invalid value '{override}' for override.")

        overriden = False
        if override != SILENT:
            for signature, _ in self.methods:
                if all(issubclass(a, b) for a, b in zip(types, signature)):
                    overriden = True
        if overriden and override == FAIL:
            raise TypeError
        elif overriden and override == WARN:
            warn('Definition ({}) overrides prior definition ({}).'.format(_fmt_t(types),
                                                                             _fmt_t(signature)),
                 TypeWarning, stacklevel=3)

        self.methods.append((types, fun))

    def add_dec(self, *types, **kwargs):
        """
        Return a decorator that adds the function it receives to the
        multimethod with the types passed as \\*args.

        You can pass in keyword argument ``override`` to control the
        overriding behavior.
        """
        self.cache = {}

        def _dec(fun):
            self.add(fun, types, kwargs.get('override', SILENT))
            return fun
        return _dec

    def __call__(self, *args, **kwargs):
        objs = self.get(*args, **kwargs)

        types = tuple(map(type, objs))

        # This code is duplicate for performance reasons.
        cached = self.cache.get(types, None)
        if cached is not None:
            return cached(*args, **kwargs)

        for signature, fun in reversed(self.methods):
            if all(issubclass(ty, sig) for ty, sig in zip(types, signature)):
                self.cache[types] = fun
                return fun(*args, **kwargs)
        raise TypeError(f'{types!r}')

    # XXX: Other Python implementations.
    def super(self, *args, **kwargs):
        """
        Like ``__call__``, only that when you give it ``super(cls, obj)``
        items, it will skip the multimethod for ``cls`` and use the one for its
        parent class.

        The normal ``__call__`` does not consider this for performance
        reasons.
        """
        objs = self.get(*args, **kwargs)
        types = tuple(
            [
                x.__thisclass__.__mro__[1] if isinstance(x, super) else type(x)
                for x in objs
            ]
        )
        nargs = [
            x.__self__ if isinstance(x, super) else x
            for x in args
        ]

        for k, elem in kwargs.items():
            if isinstance(elem, super):
                kwargs[k] = elem.__self__

        # This code is duplicate for performance reasons.
        cached = self.cache.get(types, None)
        if cached is not None:
            return cached(*nargs, **kwargs)

        for signature, fun in reversed(self.methods):
            if all(issubclass(ty, sig) for ty, sig in zip(types, signature)):
                self.cache[types] = fun
                return fun(*nargs, **kwargs)
        raise TypeError
