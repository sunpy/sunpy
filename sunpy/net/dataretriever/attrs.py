from __future__ import absolute_import

from sunpy.net.attr import Attr

class _DataSimpleAttr(Attr):
    """ A _SimpleAttr is an attribute that is not composite, i.e. that only
    has a single value, such as, e.g., Instrument('eit'). """
    def __init__(self, value):
        Attr.__init__(self)

        self.value = value

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<{cname!s}({val!r})>".format(
            cname=self.__class__.__name__, val=self.value)


class Species(_DataSimpleAttr):
    pass