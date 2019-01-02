# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

import os
from urllib.parse import urlsplit

from astropy.time import TimeDelta
from astropy.time import Time
import astropy.units as u


from datetime import timedelta
from sunpy.time import parse_time, TimeRange
from ..client import GenericClient
from sunpy import config

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['XRSClient']


class XRSClient(GenericClient):
    def _get_goes_sat_num(self, date):
        """
        Determines the satellite number for a given date.

        Parameters
        ----------

        date : `astropy.time.Time`
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
            15: TimeRange('2010-09-01', parse_time('now'))
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

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][4:]

            # 1999-01-15 as an integer.
            if int(datestamp) <= 990115:
                start = Time.strptime(datestamp, "%y%m%d")
            else:
                start = Time.strptime(datestamp, "%Y%m%d")

            almost_day = TimeDelta(1*u.day - 1*u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
        satellitenumber : int
            GOES satellite number (default = 15)
        data_type : str
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s)
        """
        # find out which satellite and datatype to query from the query times
        base_url = 'https://umbra.nascom.nasa.gov/goes/fits/'
        start_time = Time(timerange.start.strftime('%Y-%m-%d'))
        # make sure we are counting a day even if only a part of it is in the query range.
        day_range = TimeRange(timerange.start.strftime('%Y-%m-%d'),
                              timerange.end.strftime('%Y-%m-%d'))
        total_days = int(day_range.days.value) + 1
        result = list()

        # Iterate over each day in the input timerange and generate a URL for
        # it.
        for day in range(total_days):
            # It is okay to convert to datetime here as the start_time is a date
            # hence we don't necesserily gain anything.
            # This is necessary because when adding a day to a Time, we may
            # end up with the same day if the day is a leap second day
            date = start_time.datetime + timedelta(days=day)
            regex = date.strftime('%Y') + "/go{sat:02d}"
            if (date < parse_time('1999/01/15')):
                regex += date.strftime('%y%m%d') + '.fits'
            else:
                regex += date.strftime('%Y%m%d') + '.fits'
            satellitenumber = kwargs.get('satellitenumber', self._get_goes_sat_num(date))
            url = base_url + regex.format(sat=satellitenumber)
            result.append(url)
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
        chkattr = ['Time', 'Instrument', 'SatelliteNumber']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() in ('xrs', 'goes'):
                return all(chklist)
        return False
