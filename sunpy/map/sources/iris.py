# -*- coding: utf-8 -*-

from sunpy.map import GenericMap

__all__ = ['IRISMap']

class IRISMap(GenericMap):
    """
    A 2D IRIS Map
    """
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
    
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        tele = header.get('TELESCOP', '').startswith('IRIS')
        obs = header.get('INSTRUME', '').startswith('SJI')
        return tele and obs
