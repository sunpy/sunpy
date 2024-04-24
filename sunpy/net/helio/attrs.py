from sunpy.net.attr import SimpleAttr

__all__ = ['MaxRecords', 'TableName']


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class MaxRecords(SimpleAttr):
    """
    The maximum number of desired records.
    """
    def __init__(self, value):
        super().__init__(value)
        if self.value > 20000:
            raise ValueError("Helio will only return a max of 20000 results.")


class TableName(SimpleAttr):
    """
    The table to query from
    """
