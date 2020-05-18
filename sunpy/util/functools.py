"""
This file defines wrappers and variants of things in the functools standard lib.
"""

import functools

__all__ = ['seconddispatch']


def seconddispatch(func):
    """
    A variant of `functools.singledispatch` which dispatches on type of the second argument.
    """

    dispatcher = functools.singledispatch(func)

    def wrapper(*args, **kwargs):
        return dispatcher.dispatch(args[1].__class__)(*args, **kwargs)

    wrapper.dispatch = dispatcher.dispatch
    wrapper.register = dispatcher.register
    wrapper.registry = dispatcher.registry
    wrapper._clear_cache = dispatcher._clear_cache
    functools.update_wrapper(wrapper, func)

    return wrapper


# Extend the docstring with the original docs
seconddispatch.__doc__ += functools.singledispatch.__doc__
