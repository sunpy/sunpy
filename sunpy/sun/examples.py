"""
Measurement Examples
--------------------
This module provides example data of various solar related phenomena. These examples
could be used as canonical input into models. These data are taken from a variety of 
sources. 

Source
------
Sources listed in each individual example.

.. todo:: Need better sources as well as error values.

.. todo:: Create a cheat sheet function which prints out key solar values.

"""

from numpy import genfromtxt
import os
import sunpy
from sunpy.spectra import spectrum as sp

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "examples") 

def flare_spectrum():
    """Load a flare spectrum from July 23rd, 2002 GOES X5 Class flare. This is a composite
    spectrum from x-ray to gamma-rays. Spectrum below 250 keV is the 2-minute average 
    over the flare peak. Spectrum above 250 keV is the flare-integrated spectrum scaled 
    for continuity.

    Source
    ------
    RHESSI
    """    
    
    filename = os.path.abspath(os.path.join(rootdir, "flare_spectrum.txt"))
    data = genfromtxt(filename,delimiter=",", comments=';')
    
    xunit = 'Energy [keV]'
    yunit = r'Flux [photons $cm^{-2}$ s$^{-1}$ keV$^{-1}$]'
    
    # now load into spectrum object
    #spec = sp.Spectrum(data[:,0], data[:,4], width = data[:,2], name = 'RHESSI', units = [xunit, yunit])
    spec = sp.Spectrum(data[:,0], data[:,4], name = 'RHESSI', units = [xunit, yunit])
    return spec