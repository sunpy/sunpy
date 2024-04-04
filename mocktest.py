
from unittest.mock import patch
from urllib.error import HTTPError
from sunpy.net import Scraper 
from sunpy.time import TimeRange
import unittest
import pytest
import logging
log = logging.getLogger('sunpy')
log.setLevel('DEBUG')

def test_scraper_error(endpoint):
    def patch_range(self , range):
      return [f'http://httpbin.org/status/{endpoint}']
    time = TimeRange('2012/3/4', '2012/3/4 02:00') 
    pattern = "http://proba2.oma.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
    scraper = Scraper(pattern)    
    with patch.object(Scraper, 'range', new=patch_range):
       try:
        res = scraper._httpfilelist(time)
       except HTTPError as e:
         assert(e.code == endpoint)
with log.log_to_list() as log_list:  
 def test_enqueue_limit(endpoint):
  def patch_range(self , range):
    return [f'http://httpbin.org/status/{endpoint}']
  time = TimeRange('2012/3/4', '2012/3/4 02:00') 
  pattern = "http://proba2.oma.be/lyra/data/bsd/{{year:4d}}/{{month:2d}}/{{day:2d}}/{{}}_lev{{Level:1d}}_std.fits"
  scraper = Scraper(pattern)    
  with patch.object(Scraper, 'range', new=patch_range):
    try:
      res = scraper._httpfilelist(time)
    except HTTPError as e:
      assert(e.code == endpoint )
test_scraper_error(404)
test_enqueue_limit(504)



    

# Print or process the logs as needed:


