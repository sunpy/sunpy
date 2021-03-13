"""
========
Scraping
========

How to scrape data from a given `source`.

The example uses `sunpy.util.scraper.Scraper` to scrape data on a given time
range or at a specific time.
"""

from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper

#######################################################
# Include the sunpy scraper
# Include TimeRange which will be used later

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
print("Scraper Pattern: ", scraper.pattern)
print("URL now: ", scraper.now)

#######################################################
# Declare Scraper object with instrument set as 'bbso'
# and print what url does this have if try to scrape
# data for Date-time = now!

timerange = TimeRange('2006-12-01', '2006-12-01T16:00:00')
print("List of url within given time range: ", scraper.filelist(timerange))

#######################################################
# Scraping valid url for given range of date-time.

start_year = 2006
end_year = 2007
timerange = TimeRange(f'{start_year}-01-01', f'{end_year}-12-31 23:59:59.999')
urls = dict()
scraper = Scraper(url_pattern, instrument='bbso')
tarballs = scraper.filelist(timerange)

for tb_url in tarballs:
    date = scraper._extractDateURL(tb_url)
    year = date.to_datetime().year
    urls[year] = tb_url
    # print('Tarball found for year %d', year)
    # print('URL %s', tb_url)

print("Year wise arranged tarballs: ")
for each_year, tar_list in urls.items():
    print(each_year, tar_list)

#######################################################
# Search for tarballs for an year range. First get the
# list of valid url from given time range, then add it
# at the appropriate index of dictionary.

prefix = r'https://solarmonitor.org/data/'
url_pattern = prefix + r'%Y/%m/%d/fits/(\D){4}/(\D){4}_halph_fd_%Y%m%d_%H%M%S.fts.gz'

scraper = Scraper(url_pattern, regex=True)

print(scraper._URL_followsPattern("https://solarmonitor.org/data/"
      "2006/12/01/fits/bbso/bbso_halph_fd_20061201_115944.fts.gz"))

#######################################################
# Using regex for URL pattern
