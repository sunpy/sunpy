# -*- coding: utf-8 -*-
import os
import glob

from ...types import SunPyType


__all__ = []


SCHEMA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           '..', '..',
                                           'schemas',
                                           'sunpy.org',
                                           'sunpy'))


def _get_frames():
    """
    By reading the schema files, get the list of all the frames we can
    save/load.
    """
    search = os.path.join(SCHEMA_PATH, 'coordinates', 'frames', '*.yaml')
    files = glob.glob(search)

    names = []
    for fpath in files:
        path, fname = os.path.split(fpath)
        frame, _ = fname.split('-')
        exclude_schemas = []
        if frame not in exclude_schemas:
            names.append(frame)

    return names


try:
    from astropy.io.misc.asdf.tags.coordinates.frames import BaseCoordType

    __all__ = ['SunPyCoordType']

    class SunPyCoordType(BaseCoordType, SunPyType):
        _tag_prefix = "coordinates/frames/"
        name = ["coordinates/frames/" + f for f in _get_frames()]
        types = [
            'sunpy.coordinates.HeliographicCarrington',
            'sunpy.coordinates.HeliographicStonyhurst',
            'sunpy.coordinates.Heliocentric',
            'sunpy.coordinates.Helioprojective',
        ]
        requires = ['sunpy', 'astropy>=3.1']
        version = "1.0.0"

        @classmethod
        def assert_equal(cls, old, new):
            assert isinstance(new, type(old))

except ImportError:
    pass
