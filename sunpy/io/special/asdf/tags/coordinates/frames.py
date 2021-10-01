import os
import glob

import astropy.units as u
from astropy.io.misc.asdf.tags.coordinates.frames import BaseCoordType

from sunpy.coordinates import frames
from ...types import SunPyType

# Make a list of all frames with only a single schema version
sunpy_frames = list(map(lambda name: getattr(frames, name), frames.__all__))
sunpy_frames.remove(frames.HeliographicCarrington)
sunpy_frames.remove(frames.HeliographicStonyhurst)


__all__ = ['SunPyCoordType', 'HeliographicCarringtonCoordType', 'HeliographicStonyhurstCoordType']


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
        # Exclude frames with multiple schema versions
        exclude_schemas = ['heliographic_stonyhurst', 'heliographic_carrington']
        if frame not in exclude_schemas:
            names.append(frame)

    return names


class SunPyCoordType(BaseCoordType, SunPyType):
    _tag_prefix = "coordinates/frames/"
    name = ["coordinates/frames/" + f for f in _get_frames()]
    types = sunpy_frames
    requires = ['sunpy', 'astropy>=3.1']
    version = "1.0.0"

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(new, type(old))
        if new.has_data:
            assert new.data.components == old.data.components
            for comp in new.data.components:
                assert u.allclose(getattr(new.data, comp), getattr(old.data, comp))


# HeliographicStonyhurst has multiple schema versions
class HeliographicStonyhurstCoordType(SunPyCoordType):
    name = "coordinates/frames/heliographic_stonyhurst"
    types = ['sunpy.coordinates.frames.HeliographicStonyhurst']
    version = "1.1.0"
    supported_versions = ["1.0.0", "1.1.0"]


# HeliographicCarrington has multiple schema versions
class HeliographicCarringtonCoordType(SunPyCoordType):
    name = "coordinates/frames/heliographic_carrington"
    types = ['sunpy.coordinates.frames.HeliographicCarrington']
    version = "1.2.0"
    supported_versions = ["1.0.0", "1.1.0", "1.2.0"]

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        # The 1.0.0 schema should be treated as having the observer at Earth
        if cls.version == "1.0.0":
            node['frame_attributes']['observer'] = 'earth'

        return BaseCoordType.from_tree_tagged(node, ctx)
