from sunpy.net.attr import SimpleAttr

__all__ = ['MaxRecords', 'TableName']


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class MaxRecords(SimpleAttr):
    """
    The maximum number of desired records.
    """


class TableName(SimpleAttr):
    """
    The table to query from
    """
