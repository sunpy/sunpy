#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

import datetime
import urllib2
from bs4 import BeautifulSoup

from sunpy.time import parse_time, TimeRange
from sunpy.util.scraper import Scraper

from ..client import GenericClient

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['GOESClient']


class GOESClient(GenericClient):
    def _get_goes_sat_num(self, date):
        """
        Determines the satellite number for a given date.

        Parameters
        ----------

        date : `datetime.datetime`
            The date to determine which satellite is active.
        """
        goes_operational = {
            2: TimeRange('1981-01-01', '1983-04-30'),
            5: TimeRange('1983-05-02', '1984-07-31'),
            6: TimeRange('1983-06-01', '1994-08-18'),
            7: TimeRange('1994-01-01', '1996-08-13'),
            8: TimeRange('1996-03-21', '2003-06-18'),
            9: TimeRange('1997-01-01', '1998-09-08'),
            10: TimeRange('1998-07-10', '2009-12-01'),
            11: TimeRange('2006-06-20', '2008-02-15'),
            12: TimeRange('2002-12-13', '2007-05-08'),
            13: TimeRange('2006-08-01', '2006-08-01'),
            14: TimeRange('2009-12-02', '2010-10-04'),
            15: TimeRange('2010-09-01', datetime.datetime.utcnow())
        }

        results = []
        for sat_num in goes_operational:
            if date in goes_operational[sat_num]:
                # if true then the satellite with sat_num is available
                results.append(sat_num)

        if results:
            # Return the newest satellite
            return max(results)
        else:
            # if no satellites were found then raise an exception
            raise ValueError('No operational GOES satellites on {}'.format(
                date.strftime(TIME_FORMAT)))

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
        Physobs : string
            'INTENSITY', 'PARTICLE_FLUX' and 'IRRADIANCE'
        Instrument: string
            Fixed argument, 'goes'
        """
        def download_particle_flux(time_range):
            """
            Generate URLs for electron and particle fluxes measured by GOES
            for dates not older than May 30, 2015.
            URLs with `p` or `s` in their filenames stand for primary and
            secondary satellite respectively.
            """
            oldest_date = datetime.datetime(2015, 5, 30)
            total_days = (time_range.end - time_range.start).days + 1
            all_dates = time_range.split(total_days)

            url_pattern = 'ftp://ftp.swpc.noaa.gov/pub/lists/particle/{date:%Y%m%d}_G{char}_part_5m.txt'

            result = [url_pattern.format(date=day.start, char=ch)
                      for day in all_dates
                      for ch in ('p', 's') if day.start >= oldest_date]
            return result

        def download_irradiance(time_range):
            """
            Generate URLs for GOES soft x-rays flux. The result URLs will point
            to the FITS files from NASA unless required data from the last day,
            in such case the data will be obtained from NOAA.
            """
            now = datetime.datetime.utcnow()
            real_time = ((now - time_range.end).total_seconds <= 3600 * 24)

            total_days = (time_range.end - time_range.start).days + 1
            all_dates = time_range.split(total_days)

            result = list()
            # Obtain fits files within the time range
            # find out which satellite and datatype to query from the query times
            base_url = 'http://umbra.nascom.nasa.gov/goes/fits/'
            start_time = datetime.datetime.combine(timerange.start.date(),
                                                   datetime.datetime.min.time())
            total_days = int(timerange.days.value) + 1

            result = list()
            # Iterate over each day in the input timerange and generate a URL for
            # it.
            for day in range(total_days):
                date = start_time + datetime.timedelta(days=day)
                regex = "{date:%Y}/go{sat:02d}"
                if (date < parse_time('1999/01/15')):
                    regex += "{date:%y%m%d}.fits"
                else:
                    regex += "{date:%Y%m%d}.fits"
                url = base_url + regex.format(
                    date=date, sat=self._get_goes_sat_num(date))
                result.append(url)

            if real_time:  # TODO: is the above archive providing data till yesterday, or have delay?
                url_pattern = 'ftp://ftp.swpc.noaa.gov/pub/lists/xray/{date:%Y%m%d}_G{satord}_xr_{cad}m.txt'
                g = lambda url, satord, date, cad: url.format(satord=satord, date=date, cad=cad)
                result.extend([g(url_pattern, satord, now, cad) for cad in (1, 5) for satord in ('s', 'p')])
            return result

        def download_intensity(time_range):
            """
            Generate URLs for GOES SXI imager.
            """

            # TODO: update to new Scraper version #1862
            prefix = 'http://satdat.ngdc.noaa.gov/sxi/archive/fits/{sat}/{date:%Y/%m/%d/}'
            suffix = ['SXI_%Y%m%d_%H%M%S%j_{level}_{sat}.FTS'.format(pat=pat) for level in ('AB', 'BA')]
            result = list()
            total_days = (time_range.end - time_range.start).days + 1
            all_dates = time_range.split(total_days)
            for day in all_dates:
                html = urllib2.urlopen(prefix.format(date=day.end, sat=self._get_goes_sat_num(day)))
                soup = BeautifulSoup(html)
                for link in soup.findAll("a"):
                    url = str(link.get('href'))
                    for suf in suffix:
                        crawler = Scraper(suf.format(sat=self_get_goes_sat_num(day)))
                        if crawler._URL_followsPattern(url):
                            extract = crawler._extractDateURL(url)
                            valid = datetime.datetime(day.end.year, day.end.month, day.end.day,
                                                      extract.hour, extract.minute, extract.second)
                            if time_range.start <= valid <= time_range.end:
                                result.append(url)
            return result


        physobs = kwargs['physobs']
        start_dates = {'INTENSITY': datetime.datetime(2010, 6, 4), 'IRRADIANCE':datetime.datetime(2012, 9, 6),
                       'PARTICLE_FLUX': datetime.datetime(2015, 5, 30)}
        START_DATE = start_dates[physobs]
        if timerange.end < START_DATE:
            raise ValueError('Earliest date for which data is available is {:%Y-%m-%d}'.format(START_DATE))

        download_functions = {'INTENSITY': download_intensity, 'IRRADIANCE': download_irradiance,
                              'PARTICLE_FLUX': download_particle_flux}
        result = download_functions[physobs](timerange) # TODO: how to get more than one?
        return result

    def _makeimap(self):
        """
        Helper function used to hold information about source.
        """
        self.map_['source'] = 'nasa'
        self.map_['instrument'] = 'goes'
        self.map_['physobs'] = 'irradiance'
        self.map_['provider'] = 'sdac'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : list of query objects

        Returns
        -------
        boolean
            answer as to whether client can service the query
        """
        chkattr =  ['Time', 'Instrument', 'Physobs']
        physobs = ['PARTICLE_FLUX', 'IRRADIANCE', 'INTENSITY']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'goes':
                chk_var += 1
            if x.__class__.__name__ == 'Physobs' and x.value in physobs:
                chk_var += 1

        return chk_var == 2
