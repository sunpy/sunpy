# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>
# For lack of a better name, this is named magic function for now.

from __future__ import absolute_import

import inspect

def correct_argspec(fun):
    args, varargs, keywords, defaults = inspect.getargspec(fun)
    if inspect.ismethod(fun):
        args = args[1:]
    return args, varargs, keywords, defaults 


def matches_signature(fun, a, kw):
    args, varargs, keywords, defaults = correct_argspec(fun)
    if varargs is None and len(a) > len(args):
        return False
    skw = set(kw)
    sargs = set(args[len(a):])
    
    # There mayn't be unexpected parameters unless there is a **kwargs
    # in fun's signature.
    if keywords is None and skw - sargs != set():
        return False
    rest = set(args[len(a):])  - set(kw)
    
    # If there are any arguments that weren't passed but do not have
    # defaults, the signature does not match.
    defs = set() if defaults is None else set(defaults)
    if rest > defs:
        return False
    return True


class MagicFunc(object):
    def __init__(self):
        self.funcs = []
    
    def add_dec(self, condition):
        def _dec(fun):
            self.add(fun, condition)
            return fun
        return _dec
    
    def add(self, fun, condition=None):
        if condition is not None and (
            correct_argspec(fun) != correct_argspec(condition)):
            raise ValueError(
                "Signature of condition must match signature of fun."
            )
        self.funcs.append((fun, condition))
    
    def __call__(self, *args, **kwargs):
        for fun, condition in self.funcs:
            if (matches_signature(fun, args, kwargs) and
                (condition is None or condition(*args, **kwargs))):
                return fun(*args, **kwargs)
        raise TypeError
