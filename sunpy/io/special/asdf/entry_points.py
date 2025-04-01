"""
This file contains the entry points for asdf.
"""

import importlib.resources as importlib_resources

from asdf.extension import ManifestExtension
from asdf.resource import DirectoryResourceMapping


def get_resource_mappings():
    """
    Get the resource mapping instances for myschemas
    and manifests. This method is registered with the
    asdf.resource_mappings entry point.

    Returns
    -------
    list of collections.abc.Mapping
    """
    from . import resources
    resources_root = importlib_resources.files(resources)

    return [
        DirectoryResourceMapping(
            resources_root / "schemas", "asdf://sunpy.org/sunpy/schemas/"),
        DirectoryResourceMapping(
            resources_root / "manifests", "asdf://sunpy.org/sunpy/manifests/"),
    ]


def get_extensions():
    """
    Get the list of extensions.
    """
    from sunpy.io.special.asdf.converters.frames import SUNPY_FRAME_CONVERTERS
    from sunpy.io.special.asdf.converters.generic_map import GenericMapConverter

    sunpy_converters = [GenericMapConverter()] + SUNPY_FRAME_CONVERTERS
    # order here is important as asdf will pick the first matching
    # extension, start with the newest
    _manifest_uris = [
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.2.1",
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.2.0",
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.1.2",
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.1.1",
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.1.0",
        "asdf://sunpy.org/sunpy/manifests/sunpy-1.0.0",
    ]

    return [
        ManifestExtension.from_uri(uri,
                                   converters=sunpy_converters,
                                   # Register that this is a replacement for
                                   # the old extension so old files still work.
                                   # without throwing a warning.
                                   legacy_class_names=["sunpy.io.special.asdf.extension.SunpyExtension"])
        for uri
        in _manifest_uris
    ]
