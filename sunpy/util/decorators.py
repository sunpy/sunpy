# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sundry function and class decorators."""

from __future__ import print_function


import functools
import inspect
import textwrap
import types
import warnings

from sunpy.util.exceptions import SunpyDeprecationWarning
from sunpy.extern import six

__all__ = ['deprecated']


def deprecated(since, message='', name='', alternative=''):
    """
    Used to mark a function or class as deprecated.


    Parameters
    ------------
    since : str
        The release at which this API became deprecated.  This is
        required.

    message : str, optional
        Override the default deprecation message.  The format
        specifier ``func`` may be used for the name of the function,
        and ``alternative`` may be used in the deprecation message
        to insert the name of an alternative to the deprecated
        function. ``obj_type`` may be used to insert a friendly name
        for the type of object being deprecated.

    name : str, optional
        The name of the deprecated function or class; if not provided
        the name is automatically determined from the passed in
        function or class, though this is useful in the case of
        renamed functions, where the new function is just assigned to
        the name of the deprecated function.  For example::

            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative function or class name that the user may use in
        place of the deprecated object.  The deprecation warning will
        tell the user about this alternative if provided.

    """

    method_types = (classmethod, staticmethod, types.MethodType)

    def deprecate_doc(old_doc, message):
        """
        Returns a given docstring with a deprecation message prepended
        to it.
        """
        if not old_doc:
            old_doc = ''
        old_doc = textwrap.dedent(old_doc).strip('\n')
        new_doc = (('\n.. deprecated:: {since}'
                    '\n    {message}\n\n'.format(
                    **{'since': since, 'message': message.strip()})) + old_doc)
        if not old_doc:
            # This is to prevent a spurious 'unexpected unindent' warning from
            # docutils when the original docstring was blank.
            new_doc += r'\ '
        return new_doc

    def get_function(func):
        """
        Given a function or classmethod (or other function wrapper type), get
        the function object.
        """
        if isinstance(func, method_types):
            func = func.__func__
        return func

    def deprecate_function(func, message):
        """
        Returns a wrapped function that displays an
        ``AstropyDeprecationWarning`` when it is called.
        """

        if isinstance(func, method_types):
            func_wrapper = type(func)
        else:
            func_wrapper = lambda f: f

        func = get_function(func)

        def deprecated_func(*args, **kwargs):

            category = SunpyDeprecationWarning

            warnings.warn(message, category, stacklevel=2)

            return func(*args, **kwargs)

        # If this is an extension function, we can't call
        # functools.wraps on it, but we normally don't care.
        # This crazy way to get the type of a wrapper descriptor is
        # straight out of the Python 3.3 inspect module docs.
        if type(func) is not type(str.__dict__['__add__']):  # nopep8
            deprecated_func = functools.wraps(func)(deprecated_func)

        deprecated_func.__doc__ = deprecate_doc(
            deprecated_func.__doc__, message)

        return func_wrapper(deprecated_func)

    def deprecate_class(cls, message):
        """
        Returns a wrapper class with the docstrings updated and an
        __init__ function that will raise an
        ``AstropyDeprectationWarning`` warning when called.
        """
        # Creates a new class with the same name and bases as the
        # original class, but updates the dictionary with a new
        # docstring and a wrapped __init__ method.  __module__ needs
        # to be manually copied over, since otherwise it will be set
        # to *this* module (astropy.utils.misc).

        # This approach seems to make Sphinx happy (the new class
        # looks enough like the original class), and works with
        # extension classes (which functools.wraps does not, since
        # it tries to modify the original class).

        # We need to add a custom pickler or you'll get
        #     Can't pickle <class ..>: it's not found as ...
        # errors. Picklability is required for any class that is
        # documented by Sphinx.

        members = cls.__dict__.copy()

        members.update({
            '__doc__': deprecate_doc(cls.__doc__, message),
            '__init__': deprecate_function(get_function(cls.__init__),
                                           message),
        })

        return type(cls)(cls.__name__, cls.__bases__, members)

    def deprecate(obj, message=message, name=name, alternative=alternative):
        if isinstance(obj, type):
            obj_type_name = 'class'
        elif inspect.isfunction(obj):
            obj_type_name = 'function'
        elif inspect.ismethod(obj) or isinstance(obj, method_types):
            obj_type_name = 'method'
        else:
            obj_type_name = 'object'

        if not name:
            name = get_function(obj).__name__

        altmessage = ''
        if not message or type(message) is type(deprecate):
            message = ('The {func} {obj_type} is deprecated and may '
                       'be removed in a future version.')
            if alternative:
                altmessage = '\n        Use {} instead.'.format(alternative)

        message = ((message.format(**{
            'func': name,
            'name': name,
            'alternative': alternative,
            'obj_type': obj_type_name})) +
            altmessage)

        if isinstance(obj, type):
            return deprecate_class(obj, message)
        else:
            return deprecate_function(obj, message)

    if type(message) is type(deprecate):
        return deprecate(message)

    return deprecate


class add_common_docstring(object):
    """
    A function decorator that will append and/or prepend an addendum
    to the docstring of the target function.


    Parameters
    ----------
    append : `str`, optional
        A string to append to the end of the functions docstring.

    prepend : `str`, optional
        A string to prepend to the start of the functions docstring.

    **kwargs : `dict`, optional
        A dictionary to format append and prepend strings.
    """

    def __init__(self, append=None, prepend=None, **kwargs):
        if kwargs:
            append = append.format(**kwargs)
            prepend = prepend.format(**kwargs)
        self.append = append
        self.prepend = prepend

    def __call__(self, func):
        func.__doc__ = func.__doc__ if func.__doc__ else ''
        self.append = self.append if self.append else ''
        self.prepend = self.prepend if self.prepend else ''
        if self.append and isinstance(func.__doc__, six.string_types):
            func.__doc__ += self.append
        if self.prepend and isinstance(func.__doc__, six.string_types):
            func.__doc__ = self.prepend + func.__doc__
        return func
