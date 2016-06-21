from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

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

utilize sunpy/instr/tests/test_aia.py


    for units use:
    @u.quantity_input(wave=u.m)

"""

import numpy as np
import pandas as pd
import astropy.units as u
from scipy.io import readsav

import sunpy
import sunpy.data.test as test


class response():
    """
    class description here!

    Parameters
    ----------
    wave: array
        wavelengths imaged by the instrument and centered around the channel wavelength
    geoarea: float
        the geometric area of the each channel: EUV = 83.0 cm**2, UV = 30.8 cm**2
    elecperev: float
        a conversion factor for electrons per electronic volts
    ent_filter: array
        size of the entire filter on the instrument
    fp_filter: array of floats
        size of fp filter <<<<
    contam: array of floats
        maps the fluctuation of the instrument
    cross_area: array of flats
        cross area of the instrument
    secondary: array of floats
        size of secondary mirror
    primary: array
        size of the primary mirror
    ccd: array, floats
        size of the ccd


    example:
    from sunpy.instr.aia import aia                     #Q) redundant?
    effarea = aia.response.aia_inst_genx_to_dict(path_to_file)

    Notes:
    Currently the instr directory contains the file aia.py (previously aiaprep.py - I think?)
    This branch creates instr/aia directory and will have both this program aia.py and aiaprep.py

    Feedback/ Thoughts are always welcome! -Thanks!
    """

    def __init__(self, elecperev=0.0, wave=np.array, elecperdn=np.array, geoarea=0.0, ent_filter=np.array,
                 fp_filter=np.array, primary=np.array, secondary=np.array, ccd=np.array, contam=np.array,
                 cross_area=np.array, *args, **kwargs):
        # no keywords for now. emissivity, temperature, density possible later

        # parameter lists
        self.properties = ['wave', 'effarea', 'units', 'geoarea', 'scale', 'platescale',
                           'numfilters', 'wavemin', 'wavestep', 'wavenumsteps', 'wavelog',
                           'filtersincludemesh', 'meshtrans', 'usecontam', 'contamthick',
                           'usephottoelec', 'elecperev', 'usephottodn', 'elecperdn', 'useerror',
                           'fp_filter', 'ent_filter', 'primary', 'secondary', 'ccd', 'contam',
                           'cross_area']
        # TODO: want 6 EUV channels and 2 UV channels --- needs updated
        self.wavelength_centers = ['94', '131', '171', '193', '211', '304', '335']

        self.data_dictionary = {}



        # Notes: I'm trying to define these as self.properties  as something so I can define them in the definition that reads the .genx file
        #          I saw that def __repr__(self): can be used to define self.things --- need to ask about the best way to do this
        # TODO: research *args **kwargs

        # parameters to be defined from instrument file
        # self.elecperev = 0.0
        # self.wave = np.array
        # self.elecperdn = 0.0
        # self.geoarea = 0.0
        # self.ent_filter = np.array
        # self.fp_filter = np.array
        # self.primary = np.array
        # self.secondary = np.array
        # self.ccd = np.array
        # self.contam = np.array
        # self.cross_area = np.array

    def aia_inst_genx_to_dict(self, path_to_file):
        """
        This definition reads the instrument file aia_V6_all_fullinst, which was obtained from ssw_aia_response_data inside
        ssw_aia_response_genx.tar.gz (an output file saved from SolarSoft Ware (SSW)). It will extract the instrument data and save
        it a dataframe file for easier access.

         Parameters
        ----------
        path_to_file : string, the path location that leads to ssw_aia_response_data/aia_V6_all_fullinst.

        Returns
        -------
        output: dictionary, the keys are the wavelength centers for each channel and the values are dictionaries containing
        key/value pairs of the properties and values from the aia_inst_genx_.

        Notes:
        Np.recarray store information with shape(1,0) and are quite nested.

        """

        self.data_dictionary = {}

        # access np.recarray from .genx file
        ssw_array = readsav(path_to_file)
        data = ssw_array['data']

        # obtain instrument information for each channel region
        for wavelength in self.wavelength_centers:
            wave_dictionary = {}
            # pick out instrument files to look into
            for name in data.dtype.names:
                # selects instrument files inside np.recarray, eg.(A94_THICK_FULL or A94_FULL or A94_THICk_FILE)
                if name.startswith('A') and name.endswith('THICK_FULL') and wavelength in name:
                    key = wavelength

                    # store information from those files in dataframe
                    # TODO: try to refactor this to get rid of this many four loops
                    for value in data[name]:
                        for prop in self.properties:
                            wave_dictionary[prop] = [value[prop][0]]
            self.data_dictionary[wavelength] = wave_dictionary

        return self.data_dictionary

    def aia_inst_genx_to_csv(self):
        """
        This definition writes a.csv output file showing the dictionary values from ssw_genx function. The wavelength centers are the keys.
        Returns
        -------
        output: 'channel_properties.csv: dataframe with wavelength centers for each channel is a column and listed in the column are the properties from the is returned where it has keys of the properties of the channel
        """

        # write into dataframe for each channel from inner dictionaries
        df = pd.DataFrame.from_dict(self.data_dictionary)
        # save to dataframe outfile
        df.to_csv('channel_properties.csv')
        print('saved to channel_properties.csv')


    def effective_area(self, channel='A131'):
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


            # replicating the IDL version for now.     idl: if statment usephottoelec and usephottodn-- not sure why?!

            ones = np.ones(len(wave))
            print
            type(wave), type(elecperdn), type(elecperev)
            units = 12398. / wave * elecperev / elecperdn

            eff_area = ((geoarea * ent_filter * fp_filter * primary * secondary * ccd * contam) + cross_area) * units

            return eff_area

    def instrument_response(self, eff_area, gain, flat_field):
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


            # ccd == quantum efficience ==q
            # calulate gain across wavelength centers - included in master gain of aia image metadata load gain
            # contam = fluctions of the instruement
        Parameters
        ----------
        :keyword
        input: string, area from effective area calculation

        input: may need a photon-to-DN unit conversion   # digital number = counts on detectore

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
