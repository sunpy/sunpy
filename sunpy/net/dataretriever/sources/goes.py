#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

import datetime

from sunpy.time import parse_time, TimeRange

from ..client import GenericClient

__all__ = ['GOESClient']


class GOESClient(GenericClient):

    def _get_goes_sat_num(self, start, end):
        """Parses the query time to determine which GOES satellite to use."""

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
        15: TimeRange('2010-09-01', datetime.datetime.utcnow())}

        sat_list = []
        for sat_num in goes_operational:
            if ((start > goes_operational[sat_num].start and
                 start < goes_operational[sat_num].end) and
                (end > goes_operational[sat_num].start and
                 end < goes_operational[sat_num].end)):
                # if true then the satellite with sat_num is available
                sat_list.append(sat_num)

        if not sat_list:
            # if no satellites were found then raise an exception
            raise Exception('No operational GOES satellites within time range')
        else:
            return sat_list

    def _get_url_for_timerange(self, timerange, **kwargs):
        """Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
        satellite_number : int
            GOES satellite number (default = 15)
        data_type : string
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s)
        """
        # TimeRange
        if not timerange:
            return []

        start = timerange.start
        end = timerange.end
        # find out which satellite and datatype to query from the query times
        sat_num = GOESClient._get_goes_sat_num(self, start, end)
        base_url = 'http://umbra.nascom.nasa.gov/goes/fits/'

        if start < parse_time('1999/01/15'):
            url = (base_url + "%s/go%02d%s.fits") % (start.strftime("%Y"),
                sat_num[0], start.strftime("%y%m%d"))
        else:
            url = (base_url + "%s/go%02d%s.fits") % (start.strftime("%Y"),
                sat_num[0], start.strftime("%Y%m%d"))
        return [url]

    def _makeimap(self):
        """
        Helper function used to hold information about source.
        """
        self.map_['source'] = 'nasa'
        self.map_['instrument'] = 'goes'
        self.map_['phyobs'] = 'irradiance'
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
        chkattr =  ['Time', 'Instrument']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'goes':
                return all(chklist)
        return False
