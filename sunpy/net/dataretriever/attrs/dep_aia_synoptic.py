from sunpy.net.attr import SimpleAttr


__all__ = ["AIASynopticData"]


# Define a custom __dir__ to restrict tab-completion to __all__
def __dir__():
    return __all__


class AIASynopticData(SimpleAttr):
    """
    Custom attribute to indicate the use of low-resolution synoptic AIA data.

    When this attribute is present, the client will search for synoptic AIA data
    from the following source: https://jsoc1.stanford.edu/data/aia/synoptic/

    This synoptic data includes lower resolution (1k) and involves
    additional image processing steps (e.g., downsampling, time integration).

    Characteristics of the synoptic dataset can be found here:
    https://jsoc1.stanford.edu/data/aia/synoptic/README.html

    - If AIASynopticData is True, resolution defaults to 1k, overriding any user-specified resolution.
    - Defaults to True if no value is provided.

    Parameters:
        value (bool, optional): Whether to use synoptic data (defaults to True).
    """

    def __init__(self, value=True):
        if not isinstance(value, bool):
            raise ValueError("AIASynopticData value must be True or False")
        super().__init__(value)

    def __repr__(self):
        return f"AIASynopticData({self.value})"

    def is_synoptic(self):
        """Returns True if synoptic data is selected."""
        return self.value
