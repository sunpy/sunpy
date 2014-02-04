"""
Some very beta tools for IRIS
"""

import sunpy.io
import sunpy.time
import sunpy.map
from sunpy.map import GenericMap
import numpy as np

__all__ = ['SJI_to_cube', 'IRISMap']

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

def SJI_to_cube(filename, start=0, stop=None, hdu=0):
    """
    Read a SJI file and return a MapCube
    
    .. warning::
        This function is a very early beta and is not stable. Further work is 
        on going to improve SunPy IRIS support.
    
    Parameters
    ----------
    filename: string
        File to read
    
    start: int
        Temporal axis index to create MapCube from
    
    stop: int
        Temporal index to stop MapCube at
    
    hdu: int
        Choose hdu index

    Returns
    -------
    iris_cube: sunpy.map.MapCube
        A map cube of the SJI sequence
    """
    
    hdus = sunpy.io.read_file(filename)
    #Get the time delta
    time_range = sunpy.time.TimeRange(hdus[hdu][1]['STARTOBS'], hdus[hdu][1]['ENDOBS'])
    splits = time_range.split(hdus[hdu][0].shape[0])

    if not stop:
        stop = len(splits)

    headers = [hdus[hdu][1]]*(stop-start)
    datas = hdus[hdu][0][start:stop]
    
    #Make the cube:
    iris_cube = sunpy.map.Map(zip(datas,headers),cube=True)
    #Set the date/time
    for i,m in enumerate(iris_cube):
        m.meta['DATE-OBS'] = splits[i].center().isoformat()
    
    return iris_cube
