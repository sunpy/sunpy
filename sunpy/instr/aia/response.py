__author__ = 'TDWilkinson'

"""
Contains functions useful for analysing AIA data.

Goal: Use SunPy to infer plasma properties like temperature and density in multiwavelength images taken by the AIA. Two routines are necessary to calculate the response functions (while utilyzing ChiantiPy):

Wavelength response functions: calculate the amount of flux per wavelength

Temperature response functions: calculate the sensitivity of light from the plasma per temperature

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


# general format
def get_function(ion, emissivity, temperature, density, optional = None):
    """
    Statement of usefulness of function.
    Parameters
    ----------
    variable1 : what module of import is used.
        explanation
    optional: (optional) string
        A string specifying optional vairable
        e.g. strings
    """
    pass



def area():
    """
    finds the area of the instrument
    needs to work with channels
    AIA instrument response / effective area


    input: string, file (or path to file?) with instrument information

    input: a data file giving the area of the instrument

    :return: dictionary or array
    """
    pass



def temperature_response():
    """
    calculates temperature for various features using ChiantiPy


    input: string, data file

    channel: string, indicates channels used

    :return:
    """



def emissivity():
    """

    :return:
    """



def wavelength_response():
    """
    :keyword
    input: string, data file with wavelengths

    :return: float, array describing the response per wavelength
    """



def spectrum():
    """
    generate a spectral model for various solar features to be used in modules

    input:

    :return:
    """








