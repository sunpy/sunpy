from pathlib import Path

from asdf import AsdfExtension
from asdf.util import filepath_to_url

from .tags.coordinates import *  # noqa
from .tags.map import *  # noqa
from .types import SunPyType

__all__ = ['SunpyExtension']


SUNPY_SCHEMA_URI_BASE = 'http://sunpy.org/schemas/'
SCHEMA_PATH = Path(__file__).parent / "schemas"
SUNPY_URL_MAPPING = [(SUNPY_SCHEMA_URI_BASE,
                      filepath_to_url(str(SCHEMA_PATH / "sunpy.org")) + "/{url_suffix}.yaml")]


# This extension is used to register custom types that have both tags and
# schemas defined in SunPy
class SunpyExtension(AsdfExtension):
    @property
    def types(self):
        return SunPyType._tags

    @property
    def tag_mapping(self):
        return [('tag:sunpy.org:sunpy',
                 SUNPY_SCHEMA_URI_BASE + 'sunpy{tag_suffix}')]

    @property
    def url_mapping(self):
        return SUNPY_URL_MAPPING
