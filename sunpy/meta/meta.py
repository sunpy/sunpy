from __future__ import absolute_import

import numpy as np
from astropy.units import Quantity

__all__ = ['MetaItem', 'Meta']


class Meta(object):
    """
    A Meta class.  This class exposes both a dict-like interface and a
    list-like interface to metadata.

    The header may be indexed by keyword and, like a dict, the associated value
    will be returned.

    #TODO. Add ability to handle update history.
    """

    def __init__(self, metaitems):
        """
        Construct a `Meta` from an iterable.

        Parameters
        ----------
        metacards : A list of `MetaItem` objects. The items to initialize the header with.
        """
        self._items = metaitems

    def __len__(self):
        """Returns number of meta items"""
        return len(self._cards)

    def __iter__(self):
        for card in self._items:
            yield item.keyword

    def __contains__(self, metaItem):
        """Check if contains a particular metaItem"""

    def __getitem__(self, key):
        """Search for and return a metaItem with keyword"""


    def __delitem__(self, key):
        """Delete a metaItem with keyword"""

    def __repr__(self):

    def __str__(self):
        return self.tostring()

    def __eq__(self, other):
        """
        Two Headers are equal only if they have the exact same string
        representation.
        """
        return str(self) == str(other)

    def __add__(self, meta):
        """Equivalent to concatenate"""

    @property
    def items(self):
        """
        The underlying data that make up this Meta; it can be
        looked at, but it should not be modified directly.
        """
        return self._items

    @property
    def keywords(self):
        """Return all of the keywords for all of the items"""

    @property
    def values(self):
        """Returns all of the values for all of the items"""

    @property
    def comments(self):
        """View the comments associated with each keyword."""

    @property
    def sources(self):
        """View the sources of each of the items."""

    @property
    def _modified(self):
        """
        Whether or not any of the items have been modified.
        """

    @classmethod
    def fromFITS(cls, data):
        """Creates an MetaData from a FITS Header"""

    @classmethod
    def fromDict(cls, data):
        """Create a MetaData from a dictionary."""

    def to_string(cls):
        """Returns a string representation of the header"""

    def to_FITS(self):
        """Returns a FITS header"""

    @classmethod
    def copy(self, strip=False):
    """
    Make a copy of the :class:`MetaData`.
    """

    def update(self, *args, **kwargs):
    """Update the Header with new keyword values, updating the values of
    existing keywords and appending new keywords otherwise; similar to
    `dict.update`. This creates a new metaItem and adds it, keeping the old
    one for posterity."""

    def append(self, metaItems):
        """
        Appends a new item to the end of the Header, similar to `list.append`.
        """

    def extend(self, metaItems):
        """
        Appends multiple items to the header, similar to `list.extend`."""

    def remove(self, keyword):
        """
        Removes the first item with the given keyword from the header similar
        to `list.remove` if the Header object is treated as a list of keywords.

        Parameters
        ----------
        keyword : str
            The keyword of which to remove the first instance in the header

        """
        del self[self.index(keyword)]

    def is_compliant(self):
        """Check if the header is compliant (i.e. does it contain the minimum set of
        items with keywords and by checking that each item is verified."""


class MetaItem(object):
    """A metaItem represents a single item inside of a Meta object. It cannot be changed."""
    def __init__(self, keyword, value, comment, source, **kwargs):
        self._keyword = keyword
        self._value = value
        self._comment = comment
        self._modified = False
        self._source = source
        self._createDate = datetime.datetime.now()

    def __repr__(self):
        return repr((self.keyword, self.value, self.comment))

    def __str__(self):
        return self.__repr__(self)

    @property
    def createdate(self):
        """Returns the date of creation. Needed to enable history for the header"""

    @property
    def modified(self):
        """Returns true if this has been changed."""

    @property
    def source(self):
        """Returns the source of the data"""

    @property
    def keyword(self):
        """Returns the keyword name parsed from the card image."""
        return self._keyword

    @property
    def value(self):
        """The value associated with the keyword stored in this metacard."""
        return self._value

    @property
    def comment(self):
        """Get the comment attribute from the card image if not already set."""

    def _verify(self, option='warn'):
        """Verify that the card is valid. This means check that each property is not blank."""
