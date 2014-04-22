"""Hinode XRT and SOT Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = ["Jack Ireland, Jose Ivan Campos-Rozo, David Perez-Suarez"]
__email__ = "jack.ireland@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.cm import cm


__all__ = ['XRTMap', 'SOTMap']

def _lower_list(L):
    return [item.lower() for item in L]

class XRTMap(GenericMap):
    """XRT Image Map definition
    
    References
    ----------
    For a description of XRT headers
    """
    #TODO: get a link for the XRT FITS headers
    # Add in some information about the the possible filter wheel measurements
    filter_wheel1_measurements = ["Al_med", "Al_poly", "Be_med",
                                  "Be_thin", "C_poly", "Open"]
    filter_wheel2_measurements = ["Open", "Al_mesh", "Al_thick",
                                  "Be_thick", "Gband", "Ti_poly"]
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        fw1 = header.get('EC_FW1_')
        if fw1.lower() not in _lower_list(self.filter_wheel1_measurements):
            raise ValueError('Unpexpected filter wheel 1 in header.')
        fw1 = fw1.replace("_", " ")    
            
        fw2 = header.get('EC_FW2_')
        if fw2.lower() not in _lower_list(self.filter_wheel2_measurements):
            raise ValueError('Unpexpected filter wheel 2 in header.')
        fw2 = fw2.replace("_", " ")
        
        self.meta['detector'] = "XRT"
#        self.meta['instrume'] = "XRT"
        self.meta['telescop'] = "Hinode"
        
        self._name = "{0} {1}-{2}".format(self.detector, fw1, fw2)
        self._nickname = self.detector
        
        self.cmap = cm.get_cmap(name='hinodexrt')

    def _get_mpl_normalizer(self):
        """Returns a Normalize object to be used with XRT data"""
        # byte-scaled images have most likely already been scaled
        if self.dtype == np.uint8:
            return None

        mean = self.mean()
        std = self.std()
        
        vmin = max(0, mean - 3 * std)
        vmax = min(self.max(), mean + 3 * std)
        
        return colors.Normalize(vmin, vmax)


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'
        
class SOTMap(GenericMap):
	"""SOT Image Map definition
	
	References
	----------
	For a description of SOT headers
	"""
	#TODO: get a link for the SOT FITS headers
	# Add in some information about the the possible instrument, observation type,
	# observable ion and wavelength
	
	Instruments = ['SOT/WB','SOT/NB','SOT/SP','SOT/CT']

	Waves = ['6302A', 'BFI no move', 'CN bandhead 3883',
                 'Ca II H line', 'G band 4305', 'NFI no move',
                 'TF Fe I 6302', 'TF Mg I 5172', 'TF Na I 5896', 
                 'blue cont 4504', 'green cont 5550', 'red cont 6684']

	Observation_Type = ['FG (simple)', 'FG focus scan',
                            'FG shuttered I and V', 'FG shutterless I and V',
                            'FG shutterless I and V with 0.2s intervals',
                            'FG shutterless Stokes', 'SP IQUV 4D array']
	
	def __init__(self, data, header, **kwargs):
		GenericMap.__init__(self, data, header, **kwargs)

		self.meta['detector'] = "SOT"
		self.meta['telescop'] = "Hinode"
		
		self._name = self.observatory + '/' + self.instrument
		self._nickname = self.detector

                #TODO (add other options, Now all threated as intensity. This followes Hinode SDC archive)
		# StokesQUV -> grey, Velocity -> EIS, Width -> EIS, Mag Field Azi -> IDL 5 (STD gamma II)
                #'WB' -> red
		#'NB'(0 = red); (>0 = gray), # nb has 1 stokes I, the rest quv 
                #'SP' (<=1 = red); (>1 = gray) #sp has 2 stokes I, the rest quv
		color = {'SOT/WB': 'intensity', 
			 'SOT/NB': 'intensity', # For the 1st dimmension
			 'SOT/SP': 'intensity', # For the 1st 2 dimmensions
			 }

		self.cmap = cm.get_cmap('hinodesot' + color[self.instrument])


	@classmethod
	def is_datasource_for(cls, data, header, **kwargs):
		"""Determines if header corresponds to an SOT image"""
		
		return header.get('instrume') in cls.Instruments
