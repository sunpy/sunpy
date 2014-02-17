from __future__ import absolute_import

import numpy as np

from sunpy.map import GenericMap

__all__ = ['IRISMap']

class IRISMap(GenericMap):
    """
    A 2D IRIS Map
    """
    
    def __init__(self, data, header, **kwargs):    
        GenericMap.__init__(self, data, header, **kwargs)

    def iris_rot(self, missing=0.0, interpolation='bicubic', interp_param=-0.5):
        """
        Return aligned map based on header keywords
        
        Parameters
        ----------
        missing: float
           The numerical value to fill any missing points after rotation.
           Default: 0.0
        interpolation: {'nearest' | 'bilinear' | 'spline' | 'bicubic'}
            Interpolation method to use in the transform. 
            Spline uses the 
            scipy.ndimage.interpolation.affline_transform routine.
            nearest, bilinear and bicubic all replicate the IDL rot() function.
            Default: 'bicubic'
        interp_par: Int or Float
            Optional parameter for controlling the interpolation.
            Spline interpolation requires an integer value between 1 and 5 for 
            the degree of the spline fit.
            Default: 3
            BiCubic interpolation requires a flaot value between -1 and 0.
            Default: 0.5
            Other interpolation options ingore the argument.
            
        Returns
        -------
        New rotated, rescaled, translated map
        
        Notes
        -----
        Apart from interpolation='spline' all other options use a compiled 
        C-API extension. If for some reason this is not compiled correctly this
        routine will fall back upon the scipy implementation of order = 3.
        For more infomation see:
            http://sunpy.readthedocs.org/en/latest/guide/troubleshooting.html#crotate-warning
        """
        
        cords = np.matrix([[self.meta['pc1_1'], self.meta['pc1_2']],
                           [self.meta['pc2_1'], self.meta['pc2_2']]])
        center = [self.meta['CRPIX1'], self.meta['CRPIX2']]
        
        #Return a new map
        img2 = self.rotate(rmatrix=cords, rotation_center=center, recenter=False, 
                           missing=missing, interpolation=interpolation,
                           interp_param=interp_param)
             
        # modify the header to show the fact it's been corrected
        img2.meta['pc1_1'] = 1
        img2.meta['pc1_2'] = 0
        img2.meta['pc2_1'] = 0
        img2.meta['pc2_2'] = 1
    
        return img2    

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        tele = header.get('TELESCOP', '').startswith('IRIS')
        obs = header.get('INSTRUME', '').startswith('SJI')
        return tele and obs