import os
import glob

from astropy.io.misc.asdf.tags.coordinates.frames import BaseCoordType
from astropy.tests.helper import assert_quantity_allclose

from sunpy.coordinates import frames
from ...types import SunPyType

sunpy_frames = list(map(lambda name: getattr(frames, name), frames.__all__))
# Handle HeliographicCarrington separately because it has multiple schema versions
sunpy_frames.remove(frames.HeliographicCarrington)


__all__ = ['SunPyCoordType', 'HeliographicCarringtonCoordType']


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
        # Handle HeliographicCarrington separately because it has multiple schema versions
        exclude_schemas = ['heliographic_carrington']
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
                assert_quantity_allclose(getattr(new.data, comp), getattr(old.data, comp))


# Handle HeliographicCarrington specially because it has multiple schema versions
class HeliographicCarringtonCoordType(SunPyCoordType):
    name = "coordinates/frames/heliographic_carrington"
    types = ['sunpy.coordinates.frames.HeliographicCarrington']
    version = "1.1.0"
    supported_versions = ["1.0.0", "1.1.0"]

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        # The 1.0.0 schema should be treated as having the observer at Earth
        if cls.version == "1.0.0":
            node['frame_attributes']['observer'] = 'earth'

        return BaseCoordType.from_tree_tagged(node, ctx)
