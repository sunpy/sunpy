"""
========
Scraping
========

How to scrape data from a given `source`.

The example uses `sunpy.util.scraper.Scraper` to scrape data on a given time
range or at a specific time.
"""

from sunpy.util.scraper import Scraper

#######################################################
# Include the sunpy scraper

url_pattern = ('https://solarmonitor.org/data/'
               '%Y/%m/%d/fits/{instrument}/'
               '{instrument}_halph_fd_%Y%m%d_%H%M%S.fts.gz')

#######################################################
# Declare URL Pattern.
# Here %Y/%m/%d refers to Year/Month/Date.
# Similarly, %H%M%S refers to HoursMinutesSeconds.

# {some_variable} can be used to fomat the url later.
# when declaring Scraper() object you can pass
# value of this variable dynamically.

scraper = Scraper(url_pattern, instrument='bbso')
print(scraper.pattern)
print(scraper.now)

#######################################################
# Declare Scraper object with instrument set as 'bbso'
# and print what url does this have if try to scrape
# data for Date-time = now!


from sunpy.time import TimeRange
timerange = TimeRange('2006-12-01','2006-12-01T16:00:00')
print(scraper.filelist(timerange))

#######################################################
# Scraping valid url for given range of date-time.
