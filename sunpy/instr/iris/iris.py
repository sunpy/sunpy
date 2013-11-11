"""
Some very beta tools for IRIS
"""

import sunpy.io
import sunpy.time
import sunpy.map

__all__ = ['SJI_to_cube']

def SJI_to_cube(filename, start=0, stop=None):
    """
    Read a SJI file and return a MapCube
    
    ..warning::
        This function is a very early beta and is not stable. Further work is 
        on going to improve SunPy IRIS support.
    
    Parameters
    ----------
    filename: string
        File to read
    
    start:
        Temporal axis index to create MapCube from
    
    stop:
        Temporal index to stop MapCube at
    
    Returns
    -------
    
    iris_cube: sunpy.map.MapCube
        A map cube of the SJI sequence
    """
    
    hdus = sunpy.io.read_file(filename)
    #Get the time delta
    time_range = sunpy.time.TimeRange(hdus[0][1]['STARTOBS'], hdus[0][1]['ENDOBS'])
    splits = time_range.split(hdus[0][0].shape[0])

    if not stop:
        stop = len(splits)

    headers = [hdus[0][1]]*(stop-start)
    datas = hdus[0][0][start:stop]
    
    #Make the cube:
    iris_cube = sunpy.map.Map(zip(datas,headers),cube=True)
    #Set the date/time
    for i,m in enumerate(iris_cube):
        m.meta['DATE-OBS'] = splits[i].center().isoformat()
    
    return iris_cube