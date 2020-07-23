from sunpy.net.attr import SimpleAttr

__all__ = ['MaxRecords', 'TableName']


class MaxRecords(SimpleAttr):
    """
    The maximum number of desired records.
    """


class TableName(SimpleAttr):
    """
    The table to query from
    """
