
from urllib.error import HTTPError
from sunpy.net import Scraper 
from sunpy.time import TimeRange
import logging
log = logging.getLogger('sunpy')
log.setLevel('DEBUG') 
# Set logging level (optional, adjust as needed)
  # Adjust to DEBUG, INFO, WARNING, etc.

# Add a file handler to save logs (optional)
with log.log_to_file('logs.log'): 
 def test_scraper_error(endpoint):
    time = TimeRange('2012/3/4', '2012/3/4 02:00') 
    pattern = "http://proba2.oma.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
    scraper = Scraper(pattern)
    scraper.directories = [f"http://httpbin.org/{endpoint}"]
    print(scraper._httpfilelist(time))
 def test_scraper_enqueue_limit(endpoint = 504):
    time = TimeRange('2012/3/4', '2012/3/4 02:00') 
    pattern = "http://proba2.oma.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
    scraper = Scraper(pattern)
    scraper.directories = [f"http://httpbin.org/status/{endpoint}"]
    try:
     res = scraper._httpfilelist(time)
    except HTTPError as e:
     print(e)
test_scraper_error(404)
test_scraper_error(403)
test_scraper_error(400)
test_scraper_enqueue_limit()   


   


    

# Print or process the logs as needed:


