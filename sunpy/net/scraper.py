import sunpy.net
from sunpy.util.decorators import deprecated

__all__ = ['Scraper']


@deprecated("5.0", alternative="sunpy.net.Scraper")
class Scraper(sunpy.net.Scraper):
    pass


@deprecated("5.0", alternative="sunpy.net.Scraper.get_timerange_from_exdict")
def get_timerange_from_exdict(exdict):
    return sunpy.net.Scraper.get_timerange_from_exdict(exdict)


get_timerange_from_exdict.__doc__ = sunpy.net.Scraper.get_timerange_from_exdict.__doc__
