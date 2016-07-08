from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""
AIA response functions by integrating ChiantiPy

    - analyze AIA instrument response that tells efficiency of channel using instrument properties
    - then output a wavelength response
    - use spectral model (synthetic spectra) - to get ion emissivities as a function of temp/density
    - obtain temperature response based on chianti spectral contribution functions (emissivity)

AIA ccd's do not provide spectroscopic information, so there is no way to directly apply a wavelength-dependent calibration to the data.

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
from astropy.table import Table
import astropy.units as u

import os

import sunpy
import sunpy.data.test as test
import aia_read_genx2table
import aia_read_genx_test


class Response():
    """
    class description here!

    Parameters
    ----------

    Attributes
    ----------

    Notes:

    Feedback/ Thoughts are always welcome! -Thanks!
    """

    def __init__(self, channel_list, properties, **kwargs):  # data_table,

        path_to_genx_dir = kwargs.get('path_to_genx_dir', '')
        version = kwargs.get('version', 6)

        # If no data table exists, make the data table
        # if not os.path.exists('channel_properties_' + str(version) + '.csv'):
        #     aia_read_genx2table.aia_instr_properties_to_table(path_to_genx_dir, channel_list, properties, version, save=True)

        # load in astropy table  ( quicker than reading genx files each time)
        # self.data_table = Table.read('channel_properties_' + str(version) + '.csv')

        self.data_table = aia_read_genx2table.aia_instr_properties_to_table(path_to_genx_dir, channel_list, properties, version, save=True)
        self.properties = properties
        self.channel_list = channel_list




    def get_channel_data(self, channel):
        """
        This definition returns a specific row in self.data_table containing all of the channel data.
        Individual properties of the channel can be obtained by indexing the property.
        ie:  min_wavelength = get_channel_data(94)['wavemin']

        :param channel: specify which channel of which to obtain data
        :return: dictionary
        """


        if channel in self.data_table['channel']:
            index = self.channel_list.index(channel)
            channel_data = self.data_table[index]

            # Return np.array not np.ndarray
            for n, data in enumerate(channel_data):
                if type(data) == np.ndarray:
                    array = []
                    for i in data:
                        array.append(i)
                    channel_data[n] = array

            return channel_data




    def get_wavelength_range(self, channel):
        """
        wavelength range is calculated here for plotting becomes more general after calling a specific channel

        :param channel:
        :return:
        """

        num = self.get_channel_data(channel)['wavenumsteps']
        intervals = self.get_channel_data(channel)['wavestep']
        min_wave = self.get_channel_data(channel)['wavemin']

        wave = self.get_channel_data(channel)[self.properties.index('wave')]

        if channel == 'all':
            # wavelength_range = np.arange(0, float(self.data_dictionary[0][numwavesteps])) * float(
            #     (self.data_dictionary[0][wave_intervals])) + float(self.data_dictionary[0][min_wavelength])
            pass
        else:
            wavelength_range = np.arange(0, float(num)) * num + float(min_wave)

            return wavelength_range



    def effective_area(self, channel, compare_genx=False, compare_table=False):
            """
            AIA instrument response / effective area
            - contains information about the efficiency of the telescope optics.
            - trying to reformat effective area structure, convert to DN units

            Effective Area = Geo_Area * Reflectance of Primary and Secondary Mirrors * Transmission Efficiency of filters * Quantum Efficiency of CCD * Correction for time-varying contamination

            A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)
            site github ^^^     @

            Parameters
            ----------
            wavelength: int
                the wavelength center of the channel of interest

            Variables:
            ----------
            Reflectance = The combined reflectance of the primary and secondary mirrors
            Transmission_efficiency = for the focal plane and entrence filters

            Returns
            -------
            effective_area: np.array
                if compare is True: returns the effective area loaded in from the ssw genx file.
                if compare is False: returns the calculations here from Boerner


            """

            var = self.get_channel_data(channel)

            if compare_genx:
                self.eff_area = var['effarea']
            elif compare_table:
                self.eff_area = [0.312, 1.172, 2.881, 1.188, 1.206, 0.063, 0.045, 0.0192, 0.0389][
                    self.channel_list.index(channel)]
            else:
                # variables:
                wavelength = self.get_wavelength_range(channel)

                # iterate through lists for multiplication (otherwise type error occurs)
                reflectance = [i * j for i,j in zip(var['primary'], var['secondary'])]
                transmission_efficiency = [i * j for i, j in zip(var['fp_filter'], var['ent_filter'])]

                # equation:
                eff_area = var['geoarea'] * u.cm ** 2 * reflectance * transmission_efficiency * wavelength * var[
                    'contam'] * var['ccd']

                self.eff_area = eff_area

            return self.eff_area



    def wavelength_response(self, channel, compare_genx=False, compare_table=False):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature.


        Parameters
        ----------
        :keyword

        photon_to_dn?
            may need a photon-to-DN unit conversion   # digital number = counts on detector

        output: float
            outfile of instrument response per channel

        Returns
        -------
        :return: float, array describing the response per wavelength of effective area (wavelength response)
                want units cm^2 DN phot^-1
        """

        gu = u.count / u.photon
        gain_table = {94: 2.128 * gu, 131: 1.523 * gu, 171: 1.168 * gu, 195: 1.024 * gu, 211: 0.946 * gu,
                      304: 0.658 * gu, 335: 0.596 * gu, 1600: 0.125 * gu, 1700: 0.118 * gu}

        var = self.get_channel_data(channel)


        if compare_genx:
            # gives Table 2 values for instruement response
            # index = self.channel_list.index(channel)
            # self.inst_response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][index]
            # TODO: TRY  gain = elecperphot / elecperdn
            self.inst_response = self.effective_area(channel) / var['wave'] / (var['elecperphot'] / var['elecperdn'])

        elif compare_table:
            # calculations of G from Boerner Table 12
            for key, gain in gain_table.iteritems():
                if str(key) == str(channel):
                    inst_response = self.effective_area(key) / var['wave'] / gain
                    self.inst_response = inst_response # [channel]
        else:
            # Calculations of G from Boerner formulas * work in progress*
            eu = u.electron / u.count
            ccd_gain = {94: 18.3 * eu, 131: 17.6 * eu, 171: 17.7 * eu, 193: 18.3 * eu, 211: 18.3 * eu,
                        304: 18.3 * eu, 335: 17.6 * eu, 1600: 0.125 * eu, 1700: 0.118 * eu}
            for key, gain in ccd_gain.iteritems():
                if str(key) == str(channel):
                    eff_area = self.effective_area(key)
                    # g = (12398.0 / var['wave'])  / var['elecperdn'] * gain* var['elecperev']
                    # g = (12398.0*u.eV*u.angstrom / var['wave']*u.angstrom) / 3.65 /  gain
                    g = (12398.0 * u.eV * u.angstrom / var['wave'] * u.angstrom) * (u.electron / u.eV) * gain
                    instr_response = (eff_area) * g
                    self.inst_response = instr_response #[channel]

        return self.inst_response

        # NOTES:      ^^^ get the same values as the table when using all table values: self.effective_area(compare_table=True)
        #                 get one of the same values (171) when calculating it with center wavelength - fix
        #
        #     This is the
        #     closest
        #     thing in aia_bp_parse_effarea.pro - returns wavelength response
        #     effarea = thisfullstr.effarea
        #
        #
        # if KEYWORD_SET(dn) then begin
        # evperphot = 12398. / wave
        # elecperphot = (evperphot * thisfullstr.elecperev) > 1
        # elecperdn = thisfullstr.elecperdn
        # effarea = effarea * elecperphot / elecperdn
        # units = 'cm^2 DN phot^-1'
        # thisfullstr.usephottoelec = 1
        # thisfullstr.usephottodn = 1
        # endif

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
