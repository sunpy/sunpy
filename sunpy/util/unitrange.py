"""
This module provies an object that can handle a general unitrange.
"""
import astropy.units as u

__all__ = ['UnitRange']


class UnitRange:
    """
    A class to create and handle Unit ranges.

    Parameters
    ----------
    a : start of range

    b : end of range

    c : Unit of range

    """
    def __init__(self, *, a, b, unit, format=None):
        # If a is a TimeRange object, copy attributes to new instance.
        self._u  = unit
        self._t1 = a * self.unit
        self._t2 = self._t1 + b * self.unit 

    @property
    def start(self):
        """
        Get the start of unit range.
        """
        return self._t1

    @property
    def end(self):
        """
        Get the end of unit range.
        """
        return self._t2

    @property
    def unit(self):
        """
        Get the unit
        """
        return self._u

    @property
    def dt(self):
        """
        Get the length of the unit range.
        """
        return self._t2 - self._t1

    @property
    def center(self):
        """
        Gets the center of the unit range.
        """
        return self.start + self.dt / 2

    def __eq__(self, other):

        return NotImplemented

    def __ne__(self, other):

        return NotImplemented
