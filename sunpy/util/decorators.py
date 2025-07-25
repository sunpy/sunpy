"""
This module provides sunpy specific decorators.
"""
from inspect import cleandoc
from functools import wraps
from contextlib import contextmanager

import numpy as np

from astropy.utils.decorators import deprecated as _deprecated
from astropy.utils.decorators import deprecated_renamed_argument as _deprecated_renamed_argument

from sunpy.util.exceptions import SunpyDeprecationWarning

__all__ = [
    'ACTIVE_CONTEXTS',
    'deprecated',
    'deprecated_renamed_argument',
    'cached_property_based_on',
    'check_arithmetic_compatibility',
    'sunpycontextmanager',
    'add_common_docstring'
]
_NUMPY_COPY_IF_NEEDED = False if np.__version__.startswith("1.") else None
_NOT_FOUND = object()
# Stack (i.e., LIFO) of active contexts as a list of fully qualified name strings
ACTIVE_CONTEXTS = []


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
        The release at which this API became deprecated. This is
        required.
    message : str, optional
        Override the default deprecation message. The format
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
        the name of the deprecated function. For example::

            def new_function():
                ...
            oldFunction = new_function

    alternative : str, optional
        An alternative function or class name that the user may use in
        place of the deprecated object. The deprecation warning will
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


def deprecated_renamed_argument(
    old_name,
    new_name,
    since,
    arg_in_kwargs=False,
    relax=False,
    pending=False,
    warning_type=SunpyDeprecationWarning,
    alternative="",
    message="",
):
    """
    Deprecate a _renamed_ or _removed_ function argument.

    The decorator assumes that the argument with the ``old_name`` was removed
    from the function signature and the ``new_name`` replaced it at the
    **same position** in the signature. If the ``old_name`` argument is
    given when calling the decorated function the decorator will catch it and
    issue a deprecation warning and pass it on as ``new_name`` argument.

    Parameters
    ----------
    old_name : str or sequence of str
        The old name of the argument.
    new_name : str or sequence of str or None
        The new name of the argument. Set this to `None` to remove the
        argument ``old_name`` instead of renaming it.
    since : str or number or sequence of str or number
        The release at which the old argument became deprecated.
    arg_in_kwargs : bool or sequence of bool, optional
        If the argument is not a named argument (for example it
        was meant to be consumed by ``**kwargs``) set this to
        ``True``. Otherwise the decorator will throw an Exception
        if the ``new_name`` cannot be found in the signature of
        the decorated function.
        Default is ``False``.
    relax : bool or sequence of bool, optional
        If ``False`` a ``TypeError`` is raised if both ``new_name`` and
        ``old_name`` are given. If ``True`` the value for ``new_name`` is used
        and a Warning is issued.
        Default is ``False``.
    pending : bool or sequence of bool, optional
        If ``True`` this will hide the deprecation warning and ignore the
        corresponding ``relax`` parameter value.
        Default is ``False``.
    warning_type : Warning
        Warning to be issued.
        Default is `~astropy.utils.exceptions.AstropyDeprecationWarning`.
    alternative : str, optional
        An alternative function or class name that the user may use in
        place of the deprecated object if ``new_name`` is None. The deprecation
        warning will tell the user about this alternative if provided.
    message : str, optional
        A custom warning message. If provided then ``since`` and
        ``alternative`` options will have no effect.

    Raises
    ------
    TypeError
        If the new argument name cannot be found in the function
        signature and arg_in_kwargs was False or if it is used to
        deprecate the name of the ``*args``-, ``**kwargs``-like arguments.
        At runtime such an Error is raised if both the new_name
        and old_name were specified when calling the function and
        "relax=False".

    Notes
    -----
    The decorator should be applied to a function where the **name**
    of an argument was changed but it applies the same logic.

    .. warning::
        If ``old_name`` is a list or tuple the ``new_name`` and ``since`` must
        also be a list or tuple with the same number of entries. ``relax`` and
        ``arg_in_kwarg`` can be a single bool (applied to all) or also a
        list/tuple with the same number of entries like ``new_name``, etc.
    """
    return _deprecated_renamed_argument(
        old_name=old_name, new_name=new_name, since=since, arg_in_kwargs=arg_in_kwargs,
        relax=relax, pending=pending, warning_type=warning_type, alternative=alternative,
        message=message
    )


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
                # Recompute the property
                new_val = prop(instance)
                cache[prop_key] = new_val
                # Store the new attribute value after the property is computed successfully
                cache[attr_name] = new_attr_val

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
    A decorator that keeps track of active context managers in a stack.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        func_name = f"{func.__module__}.{func.__qualname__}"
        ACTIVE_CONTEXTS.append(func_name)
        gen = func(*args, **kwargs)
        value = next(gen)
        try:
            yield value
        except Exception as e:
            gen.throw(e)
        else:
            next(gen, None)
            if (removed := ACTIVE_CONTEXTS.pop()) != func_name:
                raise RuntimeError(f"Cannot remove {func_name} from tracking stack because {removed} is last active.")
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
            func.__doc__ = cleandoc(func.__doc__)  # not necessary on Python 3.13+
            func.__doc__ = func.__doc__.format(**self.kwargs)
        return func
