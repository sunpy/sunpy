# -*- coding: utf-8 -*-

from sunpy.extern import six

from asdf.asdftypes import CustomType, ExtensionTypeMeta


__all__ = ['SunPyType']


_sunpy_tags = set()


class SunPyTypeMeta(ExtensionTypeMeta):
    """
    Keeps track of `AstropyType` subclasses that are created so that they can
    be stored automatically by astropy extensions for ASDF.
    """
    def __new__(mcls, name, bases, attrs):
        cls = super(SunPyTypeMeta, mcls).__new__(mcls, name, bases, attrs)
        _sunpy_tags.add(cls)
        return cls


@six.add_metaclass(SunPyTypeMeta)
class SunPyType(CustomType):
    """
    Base class for asdf tags defined in SunPy.
    """
    organization = 'sunpy.org'
    standard = 'sunpy'
    version = '1.0.0'
