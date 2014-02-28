"""Hinode SOT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Jose Ivan Campos-Rozo"
__email__ = "hypnus1803@gmail.com"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.cm import cm


__all__ = ['SOTMap']

#~ def _lower_list(L):
	#~ return [item.lower() for item in L]

class SOTMap(GenericMap):
	"""SOT Image Map definition
	
	References
	----------
	For a description of SOT headers
	"""
	#~ #TODO: get a link for the XRT FITS headers
	#~ # Add in some information about the the possible filter wheel measurements
	
	Intruments = ['SOT/WB','SOT/NB','SOT/SP','SOT/CT']
	Waves = ['6302A', 'BFI no move', 'CN bandhead 3883', 'Ca II H line', 'G band 4305', 'NFI no move', 'TF Fe I 6302', 'TF Mg I 5172', 'TF Na I 5896', 'blue cont 4504', 'green cont 5550', 'red cont 6684']
	Observation_Type = ['FG (simple)', 'FG focus scan', 'FG shuttered I and V', 'FG shutterless I and V', 'FG shutterless I and V with 0.2s intervals', 'FG shutterless Stokes', 'SP IQUV 4D array']
	def __init__(self, data, header, **kwargs):
		
		GenericMap.__init__(self, data, header, **kwargs)
		
		
		self.meta['detector'] = "SOT"
#		self.meta['instrume'] = "SOT"
		self.meta['telescop'] = "Hinode"
		
		#self._name = self.detecto+' '+
		self._nickname = self.detector
		
		self.cmap = cm.get_cmap(name='hinodesot')

	


	@classmethod
	def is_datasource_for(cls, data, header, **kwargs):
		"""Determines if header corresponds to an SOT image"""
		instr=header.get('instrume')
		if instr in instrument:
			return instr
		#return header.get('instrume') == 'SOT/NB'
