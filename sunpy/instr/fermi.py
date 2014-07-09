import numpy as np
import copy
import math
from sunpy.time import parse_time

def nai_detector_angles():
    #returns the dictionary of Fermi/GBM NAI detector zenith and azimuth angles, in spacecraft coordinates.
    #zenith angle is measured from +z (along the LAT boresight), azimuth is measured from +x
    #see Meegan et al. (2009) for details and detector angles.
    
    #angles listed as [azimuth, zenith]
    detectors={}
    detectors['n0']  = [45.89, 20.58]
    detectors['n1']  = [45.11, 45.31]
    detectors['n2']  = [58.44, 90.21]
    detectors['n3']  = [314.87, 45.24]
    detectors['n4']  = [303.15, 90.27]
    detectors['n5']  = [3.35, 89.79]
    detectors['n6']  = [224.93, 20.43]
    detectors['n7']  = [224.62, 46.18]
    detectors['n8']  = [236.61, 89.97]
    detectors['n9']  = [135.19, 45.55]
    detectors['n10']  = [123.73, 90.42]
    detectors['n11']  = [183.74, 90.32]

    return detectors

def nai_detector_radecs(detectors, scx, scz):

    #calculates the RA/DEC for each NaI detector given spacecraft z and x RA/DEC positions
    
    scx_vector = (np.array([np.cos(np.deg2rad(scx[0]))*np.cos(np.deg2rad(scx[1])), 
                        np.sin(np.deg2rad(scx[0]))*np.cos(np.deg2rad(scx[1])),
                        np.sin(np.deg2rad(scx[1]))]))

    scz_vector = (np.array([np.cos(np.deg2rad(scz[0]))*np.cos(np.deg2rad(scz[1])), 
                        np.sin(np.deg2rad(scz[0]))*np.cos(np.deg2rad(scz[1])),
                        np.sin(np.deg2rad(scz[1]))]))

    #for each detector, do the rotation depending on the detector zenith and azimuth angles
    detector_radecs=copy.deepcopy(detectors)
    for l,d in detectors.items():
        phi = d[0]
        theta = d[1]
    
        #rotate about spacecraft z-axis first
        vx_primed = rotate_vector(scx_vector,scz_vector, np.deg2rad(phi))

        #now find spacecraft y-axis using cross product
        vy_primed = np.cross(scz_vector,vx_primed)

        #do the second part of the rotation around vy
        vz_primed = rotate_vector(scz_vector, vy_primed, np.deg2rad(theta))

        #now we should be pointing at the new RA/DEC.
        ra = np.degrees(math.atan2(vz_primed[1],vz_primed[0]))
        if ra < 0:
            ra = ra + 360.0
        dec = np.degrees(math.asin(vz_primed[2]))

        #save the RA/DEC in a dictionary
        detector_radecs[l] = [ra,dec]

    return detector_radecs
        
def rotate_vector(vector,axis,theta):
    #the Euler-Rodrigues formula for rotating vectors
    #axis is the vector to rotate around
    #theta is the angle to rotate
    #http://en.wikipedia.org/wiki/Euler-Rodrigues_parameters#Rotation_angle_and_rotation_axis
    
    axis = axis/np.sqrt(np.dot(axis,axis))
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)

    rot_matrix = np.array([[a*a+b*b-c*c-d*d, 2*(b*c+a*d),2*(b*d-a*c)],
                     [2*(b*c-a*d), a*a+c*c-b*b-d*d, 2*(c*d+a*b)],
                     [2*(b*d+a*c), 2*(c*d-a*b), a*a+d*d-b*b-c*c]])      

    #this version seems to have the wrong handedness
    #rot_matrix=np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
              #  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
              #  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    return np.dot(rot_matrix,vector)

def get_detector_separation_angles(detector_radecs,sunpos):
    #record the separation angle in degrees between the sun and each NaI detector
    angles=copy.deepcopy(detector_radecs)
    for l,d in detector_radecs.items():
        angle=separation_angle(d,sunpos)
        angles[l]=angle
    return angles

        
def separation_angle(radec1,radec2):
    #use the law of spherical cosines to calculate the separation angle between two RA/DEC positions
    #verified as correct on 2014/07/07
    cosine_of_angle = ( np.cos(np.deg2rad(90 - radec1[1])) * np.cos(np.deg2rad(90 - radec2[1])) ) + ( np.sin(np.deg2rad(90 - radec1[1])) * np.sin(np.deg2rad(90 - radec2[1])) * np.cos(np.deg2rad(radec1[0] - radec2[0])) )

    angle=np.rad2deg(np.arccos(cosine_of_angle))

    return angle

def met_to_utc(timeinsec):
    #times for GBM are in Mission Elapsed Time (MET). The reference time for this is 2001-Jan-01 00:00.
    met_ref_time = parse_time('2001-01-01 00:00')
    offset_from_utc = (met_ref_time - parse_time('1979-01-01 00:00')).total_seconds()    
    time_in_utc=parse_time(timeinsec + offset_from_utc)

    return time_in_utc
