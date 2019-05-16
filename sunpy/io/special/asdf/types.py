from asdf.types import CustomType

__all__ = ['SunPyType']


class SunPyType(CustomType):
    """
    Base class for asdf tags defined in SunPy.
    """
    organization = 'sunpy.org'
    standard = 'sunpy'
    version = '1.0.0'

    _tags = set()

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls._tags.add(cls)
