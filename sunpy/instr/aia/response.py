__author__ = 'TDWilkinson'

"""
AIA response functions by integrating ChiantiPy

    - analyze AIA instrument response then output a wavelength response
    - use spectral model (synthetic spectra)
    - obtain temperature response based on chianti spectral contribution functions



other important variables:
area
ion emissivity

electron density
temperature array

# utilize sunpy/instr/tests/test_aia.py
"""

# tools to import listed here as reference (clean out ones not used at end):
import os.path
import datetime
import csv
import copy
import socket
from itertools import dropwhile

import numpy as np
from scipy import interpolate
from scipy.integrate import trapz, cumtrapz
import astropy.units as u
import pandas as pd

import chiantipy as chpy
import chianti.core as ch

import sunpy
import sunpy.data.test as test
from sunpy.instr.aia import aiaprep
from sunpy.net import hek
from sunpy.time import parse_time
from sunpy import config
from sunpy import lightcurve
from sunpy.util.net import check_download_file
from sunpy import sun



def get_function(ion, emissivity, temperature, density, optional = None):
    """
    General format of functions where this is the statement of usefulness of function. If these were fits images, it would search the header for relevant information.

    Parameters
    ----------
    variable1 : what module of import is used.
        explanation
    optional: (optional) string
        A string specifying optional vairable
        e.g. strings

    output: tuple
    """
    pass



def effective_area():
    """
    finds the area of the instrument
    needs to work with channels
    AIA instrument response / effective area
        A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)


    input: string, file (or path to file?) with instrument information

    input: a data file giving the area of the instrument

    :return: dictionary or array
    """
    pass



def instrument_response():
    """
        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature
        R(\lambda)=A_{eff}(\lambda,t)G(\lambda)

    :keyword
    input: string, data file with wavelengths

    :return: float, array describing the response per wavelength
    """
    pass



def wavelength_response():
    """
    Wavelength response functions: calculate the amount of flux per wavelength

    :return:
    """
    pass



def emissivity():
    """

    - Use  chianti.core.mspectrum.intensityratio to get spectral line intensities and show relative emissivity per line.

    input:


    :return:
    """
    pass



def spectrum():
    """
    generate a spectral model for various solar features to be used in modules

    - develop an emissivity spectral structure based on plasma  properties
    - want to read values of physical parameters to determine feature
    - elemental abundances:  chianti.core.Spectrum.chianti.data.Abundance
        chianti.core.continuum has places to access freebound and freefree emissions

    input:

    :return:
    """
    pass



def temperature_response():
    """
    calculates temperature for various features using ChiantiPy
     - Will need to use the spectral contribution function in Chiantipy
            G(\lambda,T)
    input: string, data file

    channel: string, indicates channels used

    :return:
    """
    pass



class response()
    """
    Will I need a class at all in this?

    """
    __init__:


    pass






