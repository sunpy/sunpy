from sunpy.net.attr import SimpleAttr

__all__ = ['MaxRecords', 'Catalogue', 'TableName']


class MaxRecords(SimpleAttr):
    """
    The maximum number of desired records.
    """


class Catalogue(SimpleAttr):
    """
    The catalogue to query from.
    """


class TableName(SimpleAttr):
    """
    The table to query from
    """
