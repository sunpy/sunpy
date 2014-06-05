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

    def iris_rot(self, missing=0.0, order=3):
        """
        Return aligned map based on header keywords
        
        Parameters
        ----------
        missing: float
           The numerical value to fill any missing points after rotation.
           Default: 0.0
        order: int 0-5
            Interpolation order to be used. When using scikit-image this parameter
            is passed into :fun:`skimage.transform.warp`.
            When using scipy it is passed into 
            :fun:`scipy.ndimage.interpolation.affine_transform` where it controls 
            the order of the spline.
            
        Returns
        -------
        New rotated, rescaled, translated map
        """
        
        cords = np.matrix([[self.meta['pc1_1'], self.meta['pc1_2']],
                           [self.meta['pc2_1'], self.meta['pc2_2']]])
        
        #Return a new map
        img2 = self.rotate(rmatrix=cords, recenter=False,
                           missing=missing, order=order)
        
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
