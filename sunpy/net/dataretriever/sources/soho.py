"""SOHO Dataretriever sources subclass definitions"""

__author__ = "Ankit Kumar"
__email__  = "ankitkmr.iitk@gmail.com"

#This module was developed under funding provided by
#Google Summer of Code 2015

import urllib2
import datetime

from astropy import units as u
from bs4 import BeautifulSoup
from ..client import GenericClient

from sunpy.time import TimeRange


__all__ = ['ERNEClient']

class ERNEClient(GenericClient):

	def _get_url_for_timerange(cls, timerange, specie = 'proton'):
		"""
		Returns list of URLS to SOHO ERNE data files corresponding to value of input timerange.
		URL Source : http://srl.utu.fi/erne_data/

		The earliest data available is from 13-Feb-1996.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('1996-02-13','2013-05-15')	

		specie:  string
			Default value - proton
			Possible values - proton, alpha

		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""

		possible_specie = ['proton', 'alpha']

		#Parameter Validations
		if timerange.start < datetime.datetime(1996,02,13):
			raise ValueError('Earliest date for which SEPT data is available is 1996-02-13')

		if specie not in possible_specie:
			raise ValueError('Possible specie values: ' + '.'.join(possible_specie))


		to_continue = False
		filelists = []

		opn = urllib2.urlopen('http://srl.utu.fi/erne_data/carrot/carrot{specie}.html'.format(specie =specie[0]))

		#Getting the contents of all <tr> tags with "align" attribute having "center" value
		soup = BeautifulSoup(opn)
		results = soup.find_all("tr",{"align":"center"})
		results_string = ''

		for result in results:
			results_string = results_string + str(result)

		# Reducing the list of contents of <tr> tags to separate elemments carrying start and end date 
		# for a carrington rotation along with the carrington rotation numnber
		results_strings = results_string.split('<tr align="center">')[2:]
		final_list = [(result[5:9] + ',' + result[19:30] + ',' + result[40:51]) for result in [result[:57] for result in results_strings]]


		# Matching start and end dates of argument timerange with start and end dates of carrington rotations given on site
		# to deduce corresponding carrington numbers covered in the argument TimeRange. deduced carrington number are used 
		# to form corresponding URL to data file 
		for i in final_list[:-2]:
			carrot 	  = i[:4] 
			rot_start = i[5:16].replace(' ','-')
			if rot_start[-2] == '-':
				rot_start = rot_start[:-2] + '0' + rot_start[-1]

			rot_end   = i[17:].replace(' ','-')
			if rot_end[-2] == '-':
				rot_end = rot_end[:-2] + '0' + rot_end[-1]

			current_rotation_time = TimeRange(rot_start,rot_end)

			if (timerange.start in current_rotation_time) and (not to_continue):
				url = 'http://srl.utu.fi/erne_data/carrot/{carrot}/cr{carrot}{specie[0]}.txt'.format(carrot = carrot, specie = specie[0])
				filelists.append(url)
				if timerange.end in current_rotation_time:
					break
				else:
					to_continue = True

			if to_continue:
				url = 'http://srl.utu.fi/erne_data/carrot/{carrot}/cr{carrot}{specie}.txt'.format(carrot = carrot, specie = specie[0])
				filelists.append(url)
				if timerange.end in current_rotation_time:
					to_continue = False
					break


		return filelists

	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'university of turku'
		self.map_['instrument'] = 'soho/erne'
		self.map_['phyobs']     = 'proton and alpha particle intensities'
		self.map_['provider']   = 'space research laboratory'

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
		chkattr =  ['Time', 'Instrument']
		chklist =  [x.__class__.__name__ in chkattr for x in query]
		for x in query:
			if x.__class__.__name__ == 'Instrument' and x.value == 'soho/erne':
				return all(chklist)
		return False
