__author__ = 'TDWilkinson'

"""
AIA response functions by integrating ChiantiPy

    - analyze AIA instrument response that tells efficiency of channel using instrument properties
    - then output a wavelength response
    - use spectral model (synthetic spectra) - to get ion emissivities as a function of temp/density
    - obtain temperature response based on chianti spectral contribution functions (emissivity)


other important variables:
effective area
ion emissivity

electron density
temperature array

# utilize sunpy/instr/tests/test_aia.py
"""

import numpy as np
import pandas as pd
from scipy.io import readsav

import sunpy
import sunpy.data.test as test


def get_function(ion, emissivity, temperature, density, optional = None):
    """
    General format of functions where this is the statement of usefulness of function. If these were fits images, it would search the header for relevant information.

    Parameters
    ----------
    variable1 : what module of import is used.
        explanation
    optional: (optional) string
        A string specifying optional variable
        e.g. strings

    Returns
    -------
    output: tuple
    """
    pass
# think: how to test

def get_data_from_aia_inst_file(path_to_file):
    """
    This definition reads the instrument file aia_V6_all_fullinst file which was obtained from ssw_aia_response_data inside
    ssw_aia_response_genx.tar.gz (an output file of SolarSoft Ware (SSW)). It will extract the instrument data and save
    it a dataframe file for easier access.

     Parameters
    ----------
    :param path_to_file: string, the path location that leads to ssw_aia_response_data/aia_V6_all_fullinst.

    Returns
    -------
    output: channel_properties.csv, a comma separated dataframe file with arrays for each property per channel

    """

    # access outer dictionary from .genx file
    ssw_dict = readsav(path_to_file)

    # access inner dictionary from .genx file
    data = ssw_dict['data']

    # define wanted properties
    indexes = ['wave', 'effarea', 'units', 'geoarea', 'scale', 'platescale', 'numfilters',  'wavemin', 'wavestep', 'wavenumsteps', 'wavelog', 'filtersincludemesh', 'meshtrans', 'usecontam', 'contamthick',  'usephottoelec','elecperev','usephottodn', 'elecperdn', 'useerror', 'fp_filter', 'ent_filter', 'primary', 'secondary', 'ccd', 'contam', 'cross_area']
    channels = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']

    # write into dataframe for each channel from inner dictionaries
    df = pd.DataFrame(index=indexes, columns=channels)

    for channel in channels:

        # pick out channel files to look into
        for name in data.dtype.names:
            if name.startswith('A') and name.endswith('THICK_FULL') and channel in name:
                key = channel

                # store information from those files in dataframe
                for value in data[name]:
                    for prop in indexes:
                        df.loc[prop] = [value[prop][0]]

    # save to dataframe outfile
    df.to_csv('channel_properties.csv', sep=',')



def load_channel_as_dict(channel = ""):
    """
    This definition reads the .csv output file from the read_genx function. It will return a dictionary for each channel

     Parameters
    ----------
    :param channel, string or int, the specific channel to get properties of from the .csv file

    Returns
    -------
    output: a dictionary of one channel is returned where it has keys of the properties of the channel
    """

    dataframe = pd.DataFrame.from_csv('channel_properties.csv', sep = ',')

    # load in data for one channel and save into a dictionary
    one_channel_dict = {}
    for key,value in dataframe.iteritems():
        if type(channel) == int:
            channel = "A" + str(channel)
        if channel == str(key):
            one_channel_dict[key] = value

    assert len(one_channel_dict) != 0
    return one_channel_dict




def effective_area(channel = 'A131'):
    """
    Finds the area of the instrument
    needs to work with channels
    AIA instrument response / effective area
        A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)

    Parameters
    ----------
    filename: string, file (or path to file?) with instrument information

    input: a data file giving the area of the instrument

    channel: tell which channel to find the area for, otherwise will just find all channels.

    Returns
    -------
    effective_area: dictionary or array

    """
    if type(channel) == type:
        channel = 'A' + str(channel)

    if channel == 'all':
        pass
        # TODO: implement for all channels at once
    else:
        # load in channel properties
        channel_dictionary = load_channel_as_dict(channel)
        # define parameters  ( could use self.parameter if making a class

        for value in channel_dictionary.itervalues():
            elecperev = value['elecperev']
            wave = value['wave']
            elecperdn = value['elecperdn']
            geoarea = value['geoarea']
            ent_filter = value["ent_filter"]
            fp_filter = value['fp_filter']
            primary = value['primary']
            secondary = value['secondary']
            ccd = value['ccd']
            contam = value['contam']
            cross_area = value['cross_area']

        # replicating the IDL version for now.     idl: if statment usephottoelec and usephottodn-- not sure why?!

        ones = np.ones(len(wave))
        units = 12398.*ones / wave * elecperev / elecperdn

        eff_area = ((geoarea * ent_filter * fp_filter * primary * secondary * ccd * contam) + cross_area ) * units


        return eff_area

def instrument_response(eff_area, gain, flat_field):
    """
    Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature
        R(\lambda)=A_{eff}(\lambda,t)G(\lambda)


    notes:
    For a given position in the image plane \mathbf{x} and wavelength channel i , the pixel values can be expressed as,

    p_i(\mathbf{x})=\int_0^{\infty}\mathrm{d}\lambda\,\eta_i(\lambda)\int_{pixel\,\mathbf{x}}\mathrm{d}\theta\,I(\lambda,\theta)

    Here, \eta_i(\lambda,t,\mathbf{x}) is the efficiency function of the i^{th} channel


    effective area A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)
    gain of the CCD-camera system, G(\lambda)=(12398/\lambda/3.65)g
    flat field function F(\mathbf{x})

       Wavelength response functions: calculate the amount of flux per wavelength
     R_i(\lambda) which is equivalent to \eta_i as expressed above.


    Parameters
    ----------
    :keyword
    input: string, area from effective area calculation

    input: may need a photon-to-DN unit conversion

    output: outfile of instrument response per channel

    Returns
    -------
    :return: float, array describing the response per wavelength of effective area (wavelength response)

    """
    pass





def emissivity():
    """

    - Use  chianti.core.mspectrum.intensityratio to get spectral line intensities and show relative emissivity per line.
        chianti.core.continuum has places to access freebound and freefree emissions
    - create a line list  to show line intensities

    input:


    :return: generates chianti object - continuum - chianti model with line list
    """
    pass



def spectrum():
    """
    generate a spectral model for various solar features to be used in modules

    - develop an emissivity (chianti) spectral structure based on plasma  properties (line intensities?/ emissivity?)
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

    outfile: giving temperature response per channel

    """
    pass










