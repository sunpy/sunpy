# -*- coding: utf-8 -*-
"""
    Provides programs to process and analyze RHESSI data.

    .. warning:: This module is still in development!

"""

from __future__ import absolute_import

import urllib
import csv
from datetime import datetime
from datetime import timedelta
from itertools import dropwhile, islice, takewhile

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates

from astropy.io import fits
from astropy import units as u
from sunpy.util.odict import OrderedDict


import sunpy
import sunpy.map
import sunpy.sun.constants

from sunpy.time import TimeRange, parse_time
import sunpy.sun.constants as sun
from sunpy.sun.sun import solar_semidiameter_angular_size
from sunpy.sun.sun import sunearth_distance

__all__ = ['get_obssumm_dbase_file', 'parse_obssumm_dbase_file', 'get_obssum_filename', 'get_obssumm_file', 'parse_obssumm_file', 'backprojection']

# Measured fixed grid parameters
grid_pitch = (4.52467, 7.85160, 13.5751, 23.5542, 40.7241, 70.5309, 122.164,
              211.609, 366.646)
grid_orientation = (3.53547, 2.75007, 3.53569, 2.74962, 3.92596, 2.35647,
                    0.786083, 0.00140674, 1.57147)

data_servers = ('http://hesperia.gsfc.nasa.gov/hessidata/',
                'http://hessi.ssl.berkeley.edu/hessidata/',
                'http://soleil.i4ds.ch/hessidata/')

lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                     'navy', 'orange', 'green')

def get_obssumm_dbase_file(time_range):
    """
    Download the RHESSI observing summary database file. This file lists the
    name of observing summary files for specific time ranges.

    Parameters
    ----------
    time_range : str, TimeRange
        A TimeRange or time range compatible string

    Returns
    -------
    value : tuple
        Return a tuple (filename, headers) where filename is the local file
        name under which the object can be found, and headers is
        whatever the info() method of the object returned by urlopen.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))

    References
    ----------
    | http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing Summary Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """

    #    http://hesperia.gsfc.nasa.gov/hessidata/dbase/hsi_obssumm_filedb_200311.txt

    _time_range = TimeRange(time_range)
    data_location = 'dbase/'

    url_root = data_servers[0] + data_location
    url = url_root + _time_range.t1.strftime("hsi_obssumm_filedb_%Y%m.txt")

    f = urllib.urlretrieve(url)

    return f

def parse_obssumm_dbase_file(filename):
    """
    Parse the RHESSI observing summary database file. This file lists the
    name of observing summary files for specific time ranges along with other
    info

    Parameters
    ----------
    filename : str
        The filename of the obssumm dbase file

    Returns
    -------
    out : dict
        Return a dict containing the parsed data in the dbase file

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> f = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))
    >>> rhessi.parse_obssumm_dbase_file(f[0])

    References
    ----------
    | http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm#Observing Summary Data

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    with open(filename, "rb") as fd:
        reader = csv.reader(fd, delimiter=' ', skipinitialspace=True)
        headerline = reader.next()
        headerline = reader.next()
        headerline = reader.next()
        headerline = reader.next()

        obssumm_filename = []
        orbit_start = []
        orbit_end = []
        start_time = []
        end_time = []
        status_flag = []
        number_of_packets = []

        for row in reader:
            obssumm_filename.append(row[0])
            orbit_start.append(int(row[1]))
            orbit_end.append(int(row[2]))
            start_time.append(datetime.strptime(row[3], '%d-%b-%y'))
            end_time.append(datetime.strptime(row[5], '%d-%b-%y'))
            status_flag.append(int(row[7]))
            number_of_packets.append(int(row[8]))

        return {
            headerline[0].lower(): obssumm_filename,
            headerline[1].lower(): orbit_start,
            headerline[2].lower(): orbit_end,
            headerline[3].lower(): start_time,
            headerline[4].lower(): end_time,
            headerline[5].lower(): status_flag,
            headerline[6].lower(): number_of_packets
        }

def get_obssum_filename(time_range):
    """
    Download the RHESSI observing summary data from one of the RHESSI
    servers, parses it, and returns the name of the obssumm file relevant for
    the time range

    Parameters
    ----------
    time_range : str, TimeRange
        A TimeRange or time range compatible string

    Returns
    -------
    out : string
        Returns the filename of the observation summary file

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_filename(('2011/04/04', '2011/04/05'))

    .. note::
        This API is currently limited to providing data from whole days only.

    """
    # need to download and inspect the dbase file to determine the filename
    # for the observing summary data
    f = get_obssumm_dbase_file(time_range)
    data_location = 'metadata/catalog/'

    result = parse_obssumm_dbase_file(f[0])
    _time_range = TimeRange(time_range)

    index_number = int(_time_range.t1.strftime('%d')) - 1

    return data_servers[0] + data_location + result.get('filename')[index_number] + 's'

def get_obssumm_file(time_range):
    """
    Download the RHESSI observing summary data from one of the RHESSI
    servers.

    Parameters
    ----------
    time_range : str, TimeRange
        A TimeRange or time range compatible string

    Returns
    -------
    out : tuple
        Return a tuple (filename, headers) where filename is the local file
        name under which the object can be found, and headers is
        whatever the info() method of the object returned by urlopen.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))

    .. note::
        This API is currently limited to providing data from whole days only.

    """

    time_range = TimeRange(time_range)
    data_location = 'metadata/catalog/'

    #TODO need to check which is the closest servers
    url_root = data_servers[0] + data_location

    url = url_root + get_obssum_filename(time_range)

    print('Downloading file: ' + url)
    f = urllib.urlretrieve(url)

    return f

def parse_obssumm_file(filename):
    """
    Parse a RHESSI observation summary file.

    Parameters
    ----------
    filename : str
        The filename of a RHESSI fits file.

    Returns
    -------
    out : dict
        Returns a dictionary.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> f = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))
    >>> data = rhessi.parse_obssumm_file(f[0])

    """

    afits = fits.open(filename)
    header = afits[0].header

    reference_time_ut = parse_time(afits[5].data.field('UT_REF')[0])
    time_interval_sec = afits[5].data.field('TIME_INTV')[0]
    # label_unit = fits[5].data.field('DIM1_UNIT')[0]
    # labels = fits[5].data.field('DIM1_IDS')
    labels = ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV',
              '50 - 100 keV', '100 - 300 keV', '300 - 800 keV', '800 - 7000 keV',
              '7000 - 20000 keV']

    lightcurve_data = np.array(afits[6].data.field('countrate'))

    dim = np.array(lightcurve_data[:,0]).size

    time_array = [reference_time_ut + timedelta(0,time_interval_sec*a) for a in np.arange(dim)]

    #TODO generate the labels for the dict automatically from labels
    data = {'time': time_array, 'data': lightcurve_data, 'labels': labels}

    return header, data

def _backproject(calibrated_event_list, detector=8, pixel_size=(1.,1.), image_dim=(64,64)):
    """
    Given a stacked calibrated event list fits file create a back
    projection image for an individual detectors. This function is used by
    backprojection.

    Parameters
    ----------
    calibrated_event_list : string
        filename of a RHESSI calibrated event list
    detector : int
        the detector number
    pixel_size : 2-tuple
        the size of the pixels in arcseconds. Default is (1,1).
    image_dim : 2-tuple
        the size of the output image in number of pixels

    Returns
    -------
    out : ndarray
        Return a backprojection image.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi

    """
    afits = fits.open(calibrated_event_list)

    #info_parameters = fits[2]
    #detector_efficiency = info_parameters.data.field('cbe_det_eff$$REL')

    afits = fits.open(calibrated_event_list)

    fits_detector_index = detector + 2
    detector_index = detector - 1
    grid_angle = np.pi/2. - grid_orientation[detector_index]
    harm_ang_pitch = grid_pitch[detector_index]/1

    phase_map_center = afits[fits_detector_index].data.field('phase_map_ctr')
    this_roll_angle = afits[fits_detector_index].data.field('roll_angle')
    modamp = afits[fits_detector_index].data.field('modamp')
    grid_transmission = afits[fits_detector_index].data.field('gridtran')
    count = afits[fits_detector_index].data.field('count')

    tempa = (np.arange(image_dim[0]*image_dim[1]) %  image_dim[0]) - (image_dim[0]-1)/2.
    tempb = tempa.reshape(image_dim[0],image_dim[1]).transpose().reshape(image_dim[0]*image_dim[1])

    pixel = np.array(zip(tempa,tempb))*pixel_size[0]
    phase_pixel = (2*np.pi/harm_ang_pitch)* ( np.outer(pixel[:,0], np.cos(this_roll_angle - grid_angle)) -
                                              np.outer(pixel[:,1], np.sin(this_roll_angle - grid_angle))) + phase_map_center
    phase_modulation = np.cos(phase_pixel)
    gridmod = modamp * grid_transmission
    probability_of_transmission = gridmod*phase_modulation + grid_transmission
    bproj_image = np.inner(probability_of_transmission, count).reshape(image_dim)

    return bproj_image


def backprojection(calibrated_event_list, pixel_size=(1.,1.) * u.arcsec, image_dim=(64,64) * u.pix):
    """
    Given a stacked calibrated event list fits file create a back
    projection image.

    .. warning:: The image is not in the right orientation!

    Parameters
    ----------
    calibrated_event_list : string
        filename of a RHESSI calibrated event list
    detector : int
        the detector number
    pixel_size : `~astropy.units.Quantity` instance
        the size of the pixels in arcseconds. Default is (1,1).
    image_dim : `~astropy.units.Quantity` instance
        the size of the output image in number of pixels

    Returns
    -------
    out : RHESSImap
        Return a backprojection map.

    Examples
    --------
    >>> import sunpy.instr.rhessi as rhessi
    >>> map = rhessi.backprojection(sunpy.RHESSI_EVENT_LIST)
    >>> map.peek()

    """
    if not isinstance(pixel_size, u.Quantity):
        raise ValueError("Must be astropy Quantity in arcseconds")
    try:
       pixel_size = pixel_size.to(u.arcsec)
    except:
       raise ValueError("'{0}' is not a valid pixel_size unit".format(pixel_size.unit))
    if not (isinstance(image_dim, u.Quantity) and image_dim.unit == 'pix'):
        raise ValueError("Must be astropy Quantity in pixels")
    calibrated_event_list = sunpy.RHESSI_EVENT_LIST
    afits = fits.open(calibrated_event_list)
    info_parameters = afits[2]
    xyoffset = info_parameters.data.field('USED_XYOFFSET')[0]
    time_range = TimeRange(info_parameters.data.field('ABSOLUTE_TIME_RANGE')[0])
    
    image = np.zeros(image_dim.value)
    
    #find out what detectors were used
    det_index_mask = afits[1].data.field('det_index_mask')[0]
    detector_list = (np.arange(9)+1) * np.array(det_index_mask)
    for detector in detector_list:
        if detector > 0:
            image = image + _backproject(calibrated_event_list, detector=detector, pixel_size=pixel_size.value
										 , image_dim=image_dim.value)
    
    dict_header = {
        "DATE-OBS": time_range.center().strftime("%Y-%m-%d %H:%M:%S"),
        "CDELT1": pixel_size[0],
        "NAXIS1": image_dim[0],
        "CRVAL1": xyoffset[0],
        "CRPIX1": image_dim[0].value/2 + 0.5, 
        "CUNIT1": "arcsec",
        "CTYPE1": "HPLN-TAN",
        "CDELT2": pixel_size[1],
        "NAXIS2": image_dim[1],
        "CRVAL2": xyoffset[1],
        "CRPIX2": image_dim[0].value/2 + 0.5,
        "CUNIT2": "arcsec",
        "CTYPE2": "HPLT-TAN",
        "HGLT_OBS": 0,
        "HGLN_OBS": 0,
        "RSUN_OBS": solar_semidiameter_angular_size(time_range.center()).value,
        "RSUN_REF": sunpy.sun.constants.radius.value,
        "DSUN_OBS": sunearth_distance(time_range.center()) * sunpy.sun.constants.au.value
    }

    header = sunpy.map.MapMeta(dict_header)
    result_map = sunpy.map.Map(image, header)

    return result_map




class RHESSIFlareList(object):
    '''An object containing the RHESSI flare list, a daily updated file of every solar event detected by RHESSI.

    Attributes
    ----------
    file_list : list
        A list of dictionaries. Each dictionary contains an event from the RHESSI flare list
    meta : string
        String information describing the flare list properties in more detail
        '''
 
    def __init__(self):     
       # '''Read in the RHESSI flare list from a remote file.'''
        base_url='http://hesperia.gsfc.nasa.gov/hessidata/dbase/hessi_flare_list.txt'
        filepath=urllib.urlretrieve(base_url)

        #define the header information for this file - info in actual file not parser-friendly
        key_list = ['Flare','Start date','Start time','Peak','End','Dur (s)','Peak c/s','Total Counts','Energy (keV)',
                    'X pos (arcsec)','Y pos (arcsec)','Radial (arcsec)','AR','Flags']

        #do some manipulation of the input file
        file = open(filepath[0],'r')
        #extract useful meta information about the flare list from the footer of the input file
        metainfo = dropwhile(lambda l: l.startswith('Notes:') == False, file)
        self.meta = '\n'.join(list(metainfo))
        
        #trim the input file to ignore header and footer information
        file.seek(0)
        file_with_trimmed_header = islice(file,7,None)
        file_with_trimmed_header_and_footer = takewhile(lambda l: l.startswith('Notes:') == False,file_with_trimmed_header)
        
        #use CSV reader on the trimmed file
        csvfile = csv.reader(file_with_trimmed_header_and_footer)
        #the flare list will be a list of dictionaries
        flare_list=[]
        
        #read the file line by line
        for line in csvfile:
            #check for any blank or empty lines and ignore them 
            if not line:
                continue
            elif line[0] == '':
                continue
                
            #if not a blank line, the line should be a flare entry
            else:
                #need to split each line into its individual columns - not playing nicely with csv.reader or csv.DictReader
                flare_info=line[0].split()
                #each flare entry is its own dictionary. Use an ordered dict for better display of event information
                flare_event_dict=OrderedDict()
            
                #for each line, map the columns to the custom keys given by key_list
                #TODO - make this more pythonic
                for i in range(0,len(key_list)):
                    #Flags is the last keyword and can have multiple entries, so dump everything remaining in the line to there
                    if key_list[i] == 'Flags':
                        flare_event_dict[key_list[i]] = flare_info[i:]
                    else:
                        flare_event_dict[key_list[i]] = flare_info[i]
                #append the flare list
                flare_list.append(flare_event_dict)

        self.flare_list = flare_list



    
