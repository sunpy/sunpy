"""
This module provides SunPy specific decorators.
"""
import types
import inspect
import textwrap
import warnings
import functools
from inspect import Parameter, signature
from functools import wraps

from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyPendingDeprecationWarning

__all__ = ['deprecated']


def get_removal_version(since):
    # Work out which version this will be removed in
    since_major, since_minor = since.split('.')[:2]
    since_lts = since_minor == '0'
    if since_lts:
        major = int(since_major)
        minor = int(since_minor) + 1
    else:
        major = int(since_major) + 1
        minor = 1
    return major, minor


def deprecated(since, message='', name='', alternative='', pending=False,
               obj_type=None):
    """
    Used to mark a function or class as deprecated.

    To mark an attribute as deprecated, use `deprecated_attribute`.

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

    pending : bool, optional
        If True, uses a SunpyPendingDeprecationWarning instead of a
        ``warning_type``.

    obj_type : str, optional
        The type of this object, if the automatically determined one
        needs to be overridden.

    warning_type : `~exceptions.Warning`
        Warning to be issued.
        Default is `~sunpy.util.exceptions.SunpyDeprecationWarning`.
    """
    major, minor = get_removal_version(since)
    removal_version = f"{major}.{minor}"
    # TODO: replace this with the astropy deprecated decorator
    return _deprecated(since, message=message, name=name, alternative=alternative, pending=pending,
                       removal_version=removal_version, obj_type=obj_type,
                       warning_type=SunpyDeprecationWarning,
                       pending_warning_type=SunpyPendingDeprecationWarning)


def _deprecated(since, message='', name='', alternative='', pending=False, removal_version=None,
                obj_type=None, warning_type=SunpyDeprecationWarning,
                pending_warning_type=SunpyPendingDeprecationWarning):
    # TODO: remove this once the removal_version kwarg has been added to the upstream
    # astropy deprecated decorator
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

    def deprecate_function(func, message, warning_type=warning_type):
        """
        Returns a wrapped function that displays ``warning_type``
        when it is called.
        """

        if isinstance(func, method_types):
            func_wrapper = type(func)
        else:
            def func_wrapper(f): return f

        func = get_function(func)

        def deprecated_func(*args, **kwargs):
            if pending:
                category = pending_warning_type
            else:
                category = warning_type

            warnings.warn(message, category, stacklevel=2)

            return func(*args, **kwargs)

        # If this is an extension function, we can't call
        # functools.wraps on it, but we normally don't care.
        # This crazy way to get the type of a wrapper descriptor is
        # straight out of the Python 3.3 inspect module docs.
        if type(func) is not type(str.__dict__['__add__']):  # NOQA
            deprecated_func = functools.wraps(func)(deprecated_func)

        deprecated_func.__doc__ = deprecate_doc(
            deprecated_func.__doc__, message)

        return func_wrapper(deprecated_func)

    def deprecate_class(cls, message, warning_type=warning_type):
        """
        Update the docstring and wrap the ``__init__`` in-place (or ``__new__``
        if the class or any of the bases overrides ``__new__``) so it will give
        a deprecation warning when an instance is created.

        This won't work for extension classes because these can't be modified
        in-place and the alternatives don't work in the general case:

        - Using a new class that looks and behaves like the original doesn't
          work because the __new__ method of extension types usually makes sure
          that it's the same class or a subclass.
        - Subclassing the class and return the subclass can lead to problems
          with pickle and will look weird in the Sphinx docs.
        """
        cls.__doc__ = deprecate_doc(cls.__doc__, message)
        if cls.__new__ is object.__new__:
            cls.__init__ = deprecate_function(get_function(cls.__init__),
                                              message, warning_type)
        else:
            cls.__new__ = deprecate_function(get_function(cls.__new__),
                                             message, warning_type)
        return cls

    def deprecate(obj, message=message, name=name, alternative=alternative,
                  pending=pending, warning_type=warning_type):
        if obj_type is None:
            if isinstance(obj, type):
                obj_type_name = 'class'
            elif inspect.isfunction(obj):
                obj_type_name = 'function'
            elif inspect.ismethod(obj) or isinstance(obj, method_types):
                obj_type_name = 'method'
            else:
                obj_type_name = 'object'
        else:
            obj_type_name = obj_type

        if not name:
            name = get_function(obj).__name__

        altmessage = ''
        if not message or type(message) is type(deprecate):
            if pending:
                message = ('The {func} {obj_type} will be deprecated in '
                           'version {deprecated_version}.')
            else:
                message = ('The {func} {obj_type} is deprecated and may '
                           'be removed in {future_version}.')
            if alternative:
                altmessage = f'\n        Use {alternative} instead.'

        if removal_version is None:
            future_version = 'a future version'
        else:
            future_version = f'version {removal_version}'

        message = ((message.format(**{
            'func': name,
            'name': name,
            'deprecated_version': since,
            'future_version': future_version,
            'alternative': alternative,
            'obj_type': obj_type_name})) +
            altmessage)

        if isinstance(obj, type):
            return deprecate_class(obj, message, warning_type)
        else:
            return deprecate_function(obj, message, warning_type)

    if type(message) is type(deprecate):
        return deprecate(message)

    return deprecate


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
        func.__doc__ = func.__doc__ if func.__doc__ else ''
        self.append = self.append if self.append else ''
        self.prepend = self.prepend if self.prepend else ''
        if self.append and isinstance(func.__doc__, str):
            func.__doc__ += self.append
        if self.prepend and isinstance(func.__doc__, str):
            func.__doc__ = self.prepend + func.__doc__
        if self.kwargs:
            func.__doc__ = func.__doc__.format(**self.kwargs)
        return func


def deprecate_positional_args_since(since, keyword_only=False):
    """
    Parameters
    ----------
    since: str
        Parameter denoting last supported version.
    """
    def deprecate_positional_args(f):
        """
        Decorator for methods that issues warnings for positional arguments
        Using the keyword-only argument syntax in pep 3102, arguments after the
        * will issue a warning when passed as a positional argument.

        Parameters
        ----------
        f: function
            Function to check arguments on.

        References
        ----------
        Taken from from `scikit-learn <https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/validation.py#L1271>`__.
        Licensed under the BSD, see "licenses/SCIKIT-LEARN.rst".
        """
        nonlocal keyword_only

        sig = signature(f)
        kwonly_args = []
        all_args = []
        keyword_only = keyword_only or tuple()

        for name, param in sig.parameters.items():
            if param.kind == Parameter.POSITIONAL_OR_KEYWORD:
                all_args.append(name)
            elif param.kind == Parameter.KEYWORD_ONLY:
                kwonly_args.append(name)

        @wraps(f)
        def inner_f(*args, **kwargs):
            extra_args = len(args) - len(all_args)
            if extra_args > 0:
                for name, arg in zip(kwonly_args[:extra_args], args[-extra_args:]):
                    if name in keyword_only:
                        raise TypeError(f"{name} must be specified as a keyword argument.")

                # ignore first 'self' argument for instance methods
                args_msg = [f'{name}={arg}'
                            for name, arg in zip(kwonly_args[:extra_args],
                                                 args[-extra_args:])]
                last_supported_version = ".".join(map(str, get_removal_version(since)))
                warnings.warn(f"Pass {', '.join(args_msg)} as keyword args. "
                              f"From version {last_supported_version} "
                              "passing these as positional arguments will result in an error.",
                              SunpyDeprecationWarning)
            kwargs.update({k: arg for k, arg in zip(sig.parameters, args)})
            return f(**kwargs)
        return inner_f
    return deprecate_positional_args


_NOT_FOUND = object()


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
    attr_name :
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
        @wraps(outer)
        def inner(instance):
            """
            Parameters
            ----------
            instance:
                Any class instance that has the property ``prop``,
                and attribute ``attr``.
            """
            cache = instance.__dict__
            prop_key = prop.__name__

            # Check if our caching method has changed ouptut
            new_attr_val = getattr(instance, attr_name)
            old_attr_val = cache.get(attr_name, _NOT_FOUND)
            if old_attr_val is _NOT_FOUND or new_attr_val != old_attr_val:
                # Store the new attribute value
                cache[attr_name] = new_attr_val
                # Recompute the property
                new_val = prop(instance)
                cache[prop_key] = new_val

            return cache[prop_key]
        return inner
    return outer
