# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>

from __future__ import absolute_import

"""Provides programs to process and analyze RHESSI data.

"""

import numpy as np
import pyfits
import sunpy

# Measured fixed grid parameters
grid_pitch = [4.52467, 7.85160, 13.5751, 23.5542, 40.7241, 70.5309, 122.164, 211.609, 366.646]
grid_orientation = [3.53547, 2.75007, 3.53569, 2.74962, 3.92596, 2.35647, 0.786083, 0.00140674, 1.57147]

def backprojection(calibrated_event_list, detector=8, pixel_size=[1.,1.], image_dim=[64,64]):
    """Given a calibrated, stacked event list fits file create a back projection image."""
    calibrated_event_list = sunpy.RHESSI_EVENT_LIST
    fits = pyfits.open(calibrated_event_list)

    # here we only care about the 1st harmonic
    # hardcode this for now
    harm_ang_pitch = 1

    # Parse the fits file
    control_parameters = fits[1]
    info_parameters = fits[2]
    det_index_mask = control_parameters.data.field('det_index_mask')
    detector_efficiency = info_parameters.data.field('cbe_det_eff$$REL')    
    
    # Now get detector specific info
    fits_detector_index = detector + 2
    detector_index = detector - 1
    
    phase_map_center = fits[fits_detector_index].data.field('phase_map_ctr')
    this_detector_efficiency = detector_efficiency[0][detector_index]
    this_livetime = fits[fits_detector_index].data.field('livetime')
    this_roll_angle = fits[fits_detector_index].data.field('roll_angle')
    
    tempa = (np.arange(image_dim[0]*image_dim[1]) %  image_dim[0]) - (image_dim[0]-1)/2.
    tempb = tempa.reshape(image_dim[0],image_dim[1]).transpose().reshape(image_dim[0]*image_dim[1])
    
    pixel = np.array(zip(tempa,tempb))*pixel_size[0]
    
    phase_pixel = (2*np.pi/harm_ang_pitch)* ( np.outer(pixel[:,0], np.cos(this_roll_angle - grid_orientation[detector_index])) - np.outer(pixel[:,1], np.sin(this_roll_angle - grid_orientation[detector_index])))
    
    return phase_pixel
