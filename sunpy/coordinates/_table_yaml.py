"""
Compatibility hooks for serializing SunPy coordinate frames in Astropy tables.
"""

from astropy.coordinates import BaseCoordinateFrame
from astropy.table import serialize as table_serialize

from . import frames, metaframes


def _get_frame_classes():
    frame_classes = {}
    for module in (frames, metaframes):
        for name in module.__all__:
            cls = getattr(module, name)
            if isinstance(cls, type) and issubclass(cls, BaseCoordinateFrame):
                key = f"{cls.__module__}.{cls.__name__}"
                frame_classes[key] = cls
    return tuple(frame_classes.values())


def _coordinate_frame_representer(tag):
    def representer(dumper, obj):
        mapping = obj.info._represent_as_dict()
        return dumper.represent_mapping(tag, mapping)

    return representer


def _coordinate_frame_constructor(cls):
    def constructor(loader, node):
        mapping = loader.construct_mapping(node)
        return cls.info._construct_from_dict(mapping)

    return constructor


def _register_construct_mixin_classes(frame_classes):
    if not hasattr(table_serialize, "__construct_mixin_classes"):
        return

    existing = table_serialize.__construct_mixin_classes
    additions = tuple(
        f"{cls.__module__}.{cls.__name__}"
        for cls in frame_classes
        if f"{cls.__module__}.{cls.__name__}" not in existing
    )
    if additions:
        table_serialize.__construct_mixin_classes = existing + additions


def _register_yaml_handlers(frame_classes):
    try:
        from astropy.io.misc.yaml import AstropyDumper, AstropyLoader
    except ImportError:
        # PyYAML is optional in Astropy, so skip registration if unavailable.
        return

    for cls in frame_classes:
        tag = f"!{cls.__module__}.{cls.__name__}"
        AstropyDumper.add_representer(cls, _coordinate_frame_representer(tag))
        AstropyLoader.add_constructor(tag, _coordinate_frame_constructor(cls))


def _register_table_yaml_support():
    frame_classes = _get_frame_classes()
    _register_construct_mixin_classes(frame_classes)
    _register_yaml_handlers(frame_classes)


_register_table_yaml_support()
