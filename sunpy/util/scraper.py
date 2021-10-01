"""
This module provides a web scraper.
"""
from sunpy.net.scraper import Scraper
from sunpy.util.decorators import deprecated

__all__ = ['Scraper']


@deprecated(since="3.1", alternative="sunpy.net.Scraper")
class Scraper(Scraper):
    pass
