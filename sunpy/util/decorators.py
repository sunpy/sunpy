"""
This module provides SunPy specific decorators.
"""
from inspect import Parameter, signature
from functools import wraps
from contextlib import contextmanager

import numpy as np

from astropy.utils.decorators import deprecated as _deprecated

from sunpy.util.exceptions import SunpyDeprecationWarning, warn_deprecated

__all__ = ['deprecated', 'sunpycontextmanager', 'ACTIVE_CONTEXTS']

_NUMPY_COPY_IF_NEEDED = False if np.__version__.startswith("1.") else None
_NOT_FOUND = object()
# Public dictionary for context (de)activation tracking
ACTIVE_CONTEXTS = {}


def deprecated(
    since,
    message="",
    name="",
    alternative="",
    obj_type=None,
    warning_type=SunpyDeprecationWarning
    ):
    """
    Used to mark a function or class as deprecated.

    To mark an attribute as deprecated, use deprecated_attribute.

    Parameters
    ----------
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
    obj_type : str, optional
        The type of this object, if the automatically determined one
        needs to be overridden.
    warning_type : Warning
        Warning to be issued.
        Default is `~.SunpyDeprecationWarning`.
    """
    return _deprecated(
        since=since, message=message, name=name, alternative=alternative, pending=False,
        obj_type=obj_type, warning_type=warning_type
    )


def deprecate_positional_args_since(func=None, *, since):
    """
    Decorator for methods that issues warnings for positional arguments.

    Using the keyword-only argument syntax in pep 3102, arguments after the
    * will issue a warning when passed as a positional argument.

    Note that when you apply this, you also have to put at * in the signature
    to create new keyword only parameters!

    Parameters
    ----------
    func : callable, default=None
        Function to check arguments on.
    since : str
        The version since when positional arguments will result in error.

    Notes
    -----
    Taken from from `scikit-learn <https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/validation.py#L40>`__.
    Licensed under the BSD, see "licenses/SCIKIT-LEARN.rst".
    """
    def _inner_deprecate_positional_args(f):
        sig = signature(f)
        kwonly_args = []
        all_args = []

        for name, param in sig.parameters.items():
            if param.kind == Parameter.POSITIONAL_OR_KEYWORD:
                all_args.append(name)
            elif param.kind == Parameter.KEYWORD_ONLY:
                kwonly_args.append(name)

        @wraps(f)
        def inner_f(*args, **kwargs):
            extra_args = len(args) - len(all_args)
            if extra_args <= 0:
                return f(*args, **kwargs)

            # extra_args > 0
            args_msg = [
                f"{name}={arg}"
                for name, arg in zip(kwonly_args[:extra_args], args[-extra_args:])
            ]
            args_msg = ", ".join(args_msg)
            warn_deprecated(
                f"Pass {args_msg} as keyword args. From version "
                f"{since} passing these as positional arguments "
                "will result in an error"
            )
            kwargs.update(zip(sig.parameters, args))
            return f(**kwargs)

        return inner_f

    if func is not None:
        return _inner_deprecate_positional_args(func)

    return _inner_deprecate_positional_args


def cached_property_based_on(attr_name):
    """
    A decorator to cache the value of a property based on the output of a
    different class attribute.

    This decorator caches the values of ``getattr(instance, method)`` and
    ``prop(instance)``. When the decorated property is accessed,
    ``getattr(instance, method)`` is called. If this returns the same as its
    cached value, the cached value of ``prop`` is returned. Otherwise both
    ``meth`` and ``prop`` are recomputed, cached, and the new value of ``prop``
    is returned.

    Parameters
    ----------
    attr_name
        The name of the attribute, on which changes are checked for. The actual
        attribute is accessed using ``getattr(attr_name, instance)``.

    Notes
    -----
    The cached value of ``meth(instance)`` is stored under the key ``meth.__name__``.
    """
    def outer(prop):
        """
        prop: the property method being decorated
        """
        @wraps(prop)
        def inner(instance):
            """
            Parameters
            ----------
            instance
                Any class instance that has the property ``prop``,
                and attribute ``attr``.
            """
            cache = instance.__dict__
            prop_key = prop.__name__

            # Check if our caching method has changed output
            new_attr_val = getattr(instance, attr_name)
            old_attr_val = cache.get(attr_name, _NOT_FOUND)
            if (old_attr_val is _NOT_FOUND or
                    new_attr_val != old_attr_val or
                    prop_key not in cache):
                # Store the new attribute value
                cache[attr_name] = new_attr_val
                # Recompute the property
                new_val = prop(instance)
                cache[prop_key] = new_val

            return cache[prop_key]
        return inner
    return outer


def check_arithmetic_compatibility(func):
    """
    A decorator to check if an arithmetic operation can
    be performed between a map instance and some other operation.
    """
    # import here to reduce import complexity of `import sunpy`
    import astropy.units as u
    from astropy.nddata import NDData

    @wraps(func)
    def inner(instance, value):
        # This is explicit because it is expected that users will try to do this. This raises
        # a different error because it is expected that this type of operation will be supported
        # in future releases.
        if isinstance(value, NDData):
            return NotImplemented
        try:
            # We want to support operations between numbers and array-like objects. This includes
            # floats, ints, lists (of the aforementioned), arrays, quantities. This test acts as
            # a proxy for these possible inputs. If it can be cast to a unitful quantity, we can
            # do arithmetic with it. Broadcasting or unit mismatches are handled later in the
            # actual operations by numpy and astropy respectively.
            _ = u.Quantity(value, copy=_NUMPY_COPY_IF_NEEDED)
        except TypeError:
            return NotImplemented
        return func(instance, value)
    return inner



def sunpycontextmanager(func):
    """
    A decorator that tracks the entry and exit of a context manager,
    setting the key's value to True on entry and False on exit.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        ACTIVE_CONTEXTS[func.__name__] = True
        gen = func(*args, **kwargs)
        value = next(gen)
        try:
            yield value
        except Exception as e:
            gen.throw(e)
        else:
            next(gen, None)
            ACTIVE_CONTEXTS[func.__name__] = False
    return contextmanager(wrapper)


class add_common_docstring:
    """
    A function decorator that will append and/or prepend an addendum to the
    docstring of the target function.

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
            append = append
            prepend = prepend
        self.append = append
        self.prepend = prepend
        self.kwargs = kwargs

    def __call__(self, func):
        func.__doc__ = func.__doc__ or ''
        self.append = self.append or ''
        self.prepend = self.prepend or ''
        if self.append and isinstance(func.__doc__, str):
            func.__doc__ += self.append
        if self.prepend and isinstance(func.__doc__, str):
            func.__doc__ = self.prepend + func.__doc__
        if self.kwargs:
            func.__doc__ = func.__doc__.format(**self.kwargs)
        return func
