from sunpy.net.attr import SimpleAttr

__all__ = ["SatelliteNumber", "VersionData"]


class SatelliteNumber(SimpleAttr):
    """
    The GOES Satellite Number
    """


class VersionData(SimpleAttr):
    """
    The version of the data. The GOES 13, 14 and
    15 data has now been preprocessed but we still
    want the availability of the old data for
    reproducibility.
    """
