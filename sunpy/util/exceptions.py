"""
This module provides errors/exceptions and warnings of general use for SunPy.

Exceptions that are specific to a given package should **not** be here,
but rather in the particular package.
"""
import warnings

__all__ = ["NoMapsInFileError",
           "SunpyWarning", "SunpyUserWarning", "SunpyDeprecationWarning",
           "SunpyPendingDeprecationWarning", "SunpyMetadataWarning", "SunpyConnectionWarning",
           "warn_user", "warn_deprecated", "warn_metadata", "warn_connection"]


class NoMapsInFileError(Exception):
    """
    An error raised when a file is opened and no maps are found.
    """


class SunpyWarning(Warning):
    """
    The base warning class from which all Sunpy warnings should inherit.

    Any warning inheriting from this class is handled by the Sunpy
    logger. This warning should not be issued in normal code. Use
    "SunpyUserWarning" instead or a specific sub-class.
    """


class SunpyUserWarning(UserWarning, SunpyWarning):
    """
    The primary warning class for Sunpy.

    Use this if you do not need a specific type of warning.
    """


class SunpyMetadataWarning(UserWarning, SunpyWarning):
    """
    Warning class for cases metadata is missing.

    This does not inherit from SunpyWarning because we want to use
    stacklevel=3 to show the user where the issue occurred in their code.
    """


class SunpyDeprecationWarning(FutureWarning, SunpyWarning):
    """
    A warning class to indicate a deprecated feature.
    """


class SunpyPendingDeprecationWarning(PendingDeprecationWarning, SunpyWarning):
    """
    A warning class to indicate a soon-to-be deprecated feature.
    """


class SunpyConnectionWarning(SunpyUserWarning):
    """
    A warning class to indicate a connection warning.

    This warning should be issued when a recoverable error occurs during a connection to a remote server, such as falling back to a mirror etc.

    This will not fail the CI (via a pytest ignore) as it is not a critical warning.
    """


def warn_metadata(msg, stacklevel=1):
    """
    Raise a `SunpyMetadataWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, SunpyMetadataWarning, stacklevel + 1)


def warn_user(msg, stacklevel=1):
    """
    Raise a `SunpyUserWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, SunpyUserWarning, stacklevel + 1)


def warn_connection(msg, stacklevel=1):
    """
    Raise a `SunpyConnectionWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, SunpyConnectionWarning, stacklevel + 1)


def warn_deprecated(msg, stacklevel=1):
    """
    Raise a `SunpyDeprecationWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, SunpyDeprecationWarning, stacklevel + 1)
