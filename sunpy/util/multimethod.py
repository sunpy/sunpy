# Copyright (c) 2011 Florian Mayer <florian.mayer@bitsrc.org>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""
Multimethod implementation in pure Python.
"""

from __future__ import absolute_import, division, print_function

from warnings import warn

from sunpy.extern.six.moves import zip, map
from sunpy.extern import six

__all__ = ['TypeWarning', 'MultiMethod']

SILENT = 0
WARN = 1
FAIL = 2

def _fmt_t(types):
    return ', '.join(type_.__name__ for type_ in types)


class TypeWarning(UserWarning):
    pass


class MultiMethod(object):
    """ A multimethod is a callable object that decides which code to execute
    based on the type of one or more of its arguments.

    Parameters
    ----------
    get : function
        function which receives args and kwargs and returns a tuple of
        values to consider for dispatch.
    """
    def __init__(self, get):
        self.get = get

        self.methods = []
        self.cache = {}

    def add(self, fun, types, override=SILENT):
        """ Add fun to the multimethod. It will be executed if get returns
        values of the types passed as types. Must return tuples of same
        length for any input.

        Parameters
        ----------
        fun : function
            function to be added to the multimethod
        types : tuple of classes
            types for which the function is executed
        override : SILENT, WARN or FAIL
            control behaviour when overriding existing definitions.
            If it is set to SILENT, prior definitions are silently
            overridden, if it is set to WARN a TypeWarning
            will be issued, and with FAIL a TypeError is raised when
            attempting to override an existing definition.
        """
        overriden = False
        if override:
            for signature, _ in self.methods:
                if all(issubclass(a, b) for a, b in zip(types, signature)):
                    overriden = True
        if overriden and override == FAIL:
            raise TypeError
        elif overriden and override == WARN:
            # pylint: disable=W0631
            warn(
                'Definition ({0}) overrides prior definition ({1}).'.format(
                _fmt_t(types), _fmt_t(signature)),
                TypeWarning,
                stacklevel=3
            )
        elif overriden:
            raise ValueError('Invalid value for override.')
        self.methods.append((types, fun))

    def add_dec(self, *types, **kwargs):
        """ Return a decorator that adds the function it receives to the
        multimethod with the types passed as \*args. Using keyword arg
        override to control overriding behaviour. Compare add.
        """
        self.cache = {}
        def _dec(fun):
            self.add(fun, types, kwargs.get('override', SILENT))
            return fun
        return _dec

    def __call__(self, *args, **kwargs):
        objs = self.get(*args, **kwargs)

        # pylint: disable=W0141
        types = tuple(map(type, objs))

        # This code is duplicate for performance reasons.
        cached = self.cache.get(types, None)
        if cached is not None:
            return cached(*args, **kwargs)

        for signature, fun in reversed(self.methods):
            if all(issubclass(ty, sig) for ty, sig in zip(types, signature)):
                self.cache[types] = fun
                return fun(*args, **kwargs)
        raise TypeError('{0!r}'.format(types))

    # XXX: Other Python implementations.
    def super(self, *args, **kwargs):
        """ Like __call__, only that when you give it super(cls, obj) items,
        it will skip the multimethod for cls and use the one for its parent
        class. The normal __call__ does not consider this for performance
        reasons. """
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

        for k, elem in six.iteritems(kwargs):
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
