#!/usr/bin/python
# -*- coding: utf8 -*-

####################### Function diff_rot  ###########################################################################
#               Written by Jose Ivan Campos Rozo                                                                     #
#               Date: August 16 2012                                                                                 #
#               Reference:                                                                                           #
#                       http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro                                 #
#               This function let us computes the differential rotation of the sun                                   # 
#               Inputs:                                                                                              #          
#                      ddays: Number of days that I want to rotate                                                   #
#                      latitude: heliographic coordinate latitude in Degrees                                         #
#                                                                                                                    #
#               Optional Inputs:                                                                                     #
#                       allen: Use values from Allen, Astrophysical Quantities                                       #
#                       howard: Use values for small magnetic features from Howard et al.                            #
#                       synodic: Use synodic rotation rate                                                           #
#                       sidereal: Use sidereal rotation rate                                                         #
#               Using the function:                                                                                  #
#                       Examples:                                                                                    #
#                               rotation=diff_rot(2,30) %without optional inputs (Default) where 2 is                #
#                                                       % number of the days and 30 is latitude in degrees.          #
#                               rotation=diff_rot(2,30,'allen')  % With optional input allen, you must               #
#                                                                % realize that the optional input must be           #
#                                                                % a string.                                         #
#                               rotation=diff_rot(2,30,'howard') % This optional input use values from Howard        #
#                                                                % and it's equal to without optional input.         #
#                               rotation=diff_rot(2,30,'sidereal')  %  Default                                       #
#                               rotation=diff_rot(2,30,'synodic') %  Use sidereal rotation rate to calculate         #
#                                                                 %  the rotation                                    #
#                                                                                                                    #
#               Output: This function computes the change in longitude over days (ddays) in Degrees                  #
#                                                                                                                    #
#                                                                                                                    #
#                                                                                                                    #
#                                                                                                                    #
######################################################################################################################
from math import pi, sin, cos
from numpy import deg2rad as d2r
import numpy as np

def diff_rot(ddays,latitude,optional=None):
    """
    This function computes the change in longitude over days in degrees,
    diff_rot(ddays,latitud,[optional])
        ddays: Number of days that I want to rotate.
        latitude: heliographic coordinate latitude in Degrees.
        Optional inputs:
        ----------------
                         Howard: Use values for small magnetic
                                 features from Howard et al.
                         synodic: Use synodic rotation rate.
                         sidereal: Use sidereal rotation rate.
                         allen: Use values from Allen, Astrophysical Quantities
    
        Returns:
        -------
                return the change in longitude over days (units=degrees)
        
    """

    sin2l = (sin(d2r(latitude)))**2
    sin4l = sin2l**2
    if (optional==None) or (optional=='howard') or (optional=='sidereal')or (optional=='synodic'):
	    rotation=1.*(10**(-6))*ddays*(2.894-0.428*sin2l-0.37*sin4l)*24.*3600./d2r(1)
	    if optional=='synodic':
		    rotation = rotation-0.9856*ddays
			 
		 
    if optional=='allen':
	    rotation= ddays*(14.44-(3.0*sin2l))

    return np.round(rotation,4)
         
    


     

     
