"""SOHO Dataretriever sources subclass definitions"""

__author__ = "Ankit Kumar"
__email__ = "ankitkmr.iitk@gmail.com"

#This module was developed under funding provided by
#Google Summer of Code 2015

import urllib2
import datetime

from bs4 import BeautifulSoup
from ..client import GenericClient

from sunpy.time import TimeRange

__all__ = ['ERNEClient']


class ERNEClient(GenericClient):
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS to SOHO ERNE data files corresponding to value of input timerange and species.
        URL Source : http://srl.utu.fi/erne_data/

        The earliest data available is from 13-Feb-1996.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
            Example value -  TimeRange('1996-02-13','2013-05-15')

        species:  string
            Default value - proton
            Possible values - proton, alpha

        Returns
        -------
        urls : list
            list of URLs corresponding to the requested time range

        Examples
        --------
        >>> from sunpy.time.timerange import TimeRange
        >>> from sunpy.net.vso.attrs import Time, Instrument
        >>> import sunpy.net.dataretriever.sources.soho as soho
        >>> LCClient = soho.ERNEClient()

        >>> qr1 = LCClient.query(Time(TimeRange('2003-03-01','2003-04-04')), Instrument('erne'), Species('alpha'))
        >>> res = LCClient.get(qr1)
        >>> download_list = res.wait()

        References
        ----------
        They are available at the srl server,
        | http://srl.utu.fi/erne_data/

        and have file names of type
        | http://srl.utu.fi/erne_data/carrot/1906/cr1906p.txt
        | http://srl.utu.fi/erne_data/carrot/1906/cr1906a.txt

        """

        #Parameter Validations
        if timerange.start < datetime.datetime(1996, 02, 13):
            raise ValueError(
                'Earliest date for which SEPT data is available is 1996-02-13')

        if 'species' in kwargs:
            species = kwargs['species']
        else:
             raise ValueError(
                'No species defined: alpha or protons')

        to_continue = False
        filelists = []

        opn = urllib2.urlopen(
            'http://srl.utu.fi/erne_data/carrot/carrot{species}.html'.format(
                species=species[0]))

        #Getting the contents of all <tr> tags with "align" attribute having "center" value
        soup = BeautifulSoup(opn)
        results = soup.find_all("tr", {"align": "center"})
        results_string = ''

        for result in results:
            results_string = results_string + str(result)

        # Reducing the list of contents of <tr> tags to separate elemments carrying start and end date
        # for a carrington rotation along with the carrington rotation numnber
        results_strings = results_string.split('<tr align="center">')[2:]
        final_list = [(result[5:9] + ',' + result[19:30] + ',' + result[40:51])
                      for result in [result[:57]
                                     for result in results_strings]]

        # Matching start and end dates of argument timerange with start and end dates of carrington rotations given on site
        # to deduce corresponding carrington numbers covered in the argument TimeRange. deduced carrington number are used
        # to form corresponding URL to data file
        for i in final_list[:-2]:
            carrot = i[:4]
            rot_start = i[5:16].replace(' ', '-')
            if rot_start[-2] == '-':
                rot_start = rot_start[:-2] + '0' + rot_start[-1]

            rot_end = i[17:].replace(' ', '-')
            if rot_end[-2] == '-':
                rot_end = rot_end[:-2] + '0' + rot_end[-1]

            current_rotation_time = TimeRange(rot_start, rot_end)

            if (timerange.start in current_rotation_time) and (
                    not to_continue):
                url = 'http://srl.utu.fi/erne_data/carrot/{carrot}/cr{carrot}{species[0]}.txt'.format(
                    carrot=carrot,
                    species=species[0])
                filelists.append(url)
                if timerange.end in current_rotation_time:
                    break
                else:
                    to_continue = True

            if to_continue:
                url = 'http://srl.utu.fi/erne_data/carrot/{carrot}/cr{carrot}{species}.txt'.format(
                    carrot=carrot,
                    species=species[0])
                filelists.append(url)
                if timerange.end in current_rotation_time:
                    to_continue = False
                    break

        return filelists

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'soho'
        self.map_['instrument'] = 'erne'
        self.map_['phyobs'] = 'proton and alpha particle intensities'
        self.map_['provider'] = 'university of turku'

    @classmethod
    def _can_handle_query(cls, *query, **kwargs):
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
        chk_var = 0
        for x in query:
            if (x.__class__.__name__ == 'Instrument' and
                x.value.lower() == 'erne'):

                chk_var += 1

            elif x.__class__.__name__ == 'Species' and x.value in ['alpha', 'proton']:
                chk_var += 1

        if(chk_var == 2):
            return True
        return False