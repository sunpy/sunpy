"""
This file contains the entry points for asdf.
"""
import sys

from asdf.extension import ManifestExtension
from asdf.resource import DirectoryResourceMapping

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources


def get_resource_mappings():
    """
    Get the resource mapping instances for myschemas
    and manifests.  This method is registered with the
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

    return [
        ManifestExtension.from_uri("asdf://sunpy.org/sunpy/manifests/sunpy-1.0.0",
                                   converters=sunpy_converters,
                                   # Register that this is a replacement for
                                   # the old extension so old files still work.
                                   # without throwing a warning.
                                   legacy_class_names=["sunpy.io.special.asdf.extension.SunpyExtension"])
    ]
