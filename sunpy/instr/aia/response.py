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

        self.data_table = aia_read_genx2table.aia_instr_properties_to_table(path_to_genx_dir, channel_list, properties,
                                                                            version, save=True)
        self.properties = properties
        self.channel_list = channel_list

    def get_channel_data(self, channel):
        """
        This definition returns a specific row in self.data_table containing all of the channel data.
        Individual properties of the channel can be obtained by indexing the property.
        ie:  min_wavelength = get_channel_data(94)['wavemin']

        :param channel: specify which channel of which to obtain data
        :return: Table row that acts like dictionary
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
            wavelength_range = np.arange(0, float(num)) * float(intervals) + float(min_wave)

            return wavelength_range

    def effective_area(self, channel, total_range=False, compare_genx=False, print_comparison_ratio=False):
        """
        AIA instrument response / effective area
        - contains information about the efficiency of the telescope optics.
        - trying to reformat effective area structure, convert to DN units

        Effective Area = Geo_Area * Reflectance of Primary and Secondary Mirrors * Transmission Efficiency of filters * Quantum Efficiency of CCD * Correction for time-varying contamination

        # formula from paper:
        A eff (wavelength, t) = Ageo*Rp(wavelength)*Rs(wavelength)*Te(wavelength)*Tf(wavelength)*D(wavelength, t)*Q(wavelength).

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


        total range returns entire wavelength range
        otherwise returns integer value for specific channel
        """

        var = self.get_channel_data(channel)

        if compare_genx:
            self.eff_area = var['effarea'] * u.cm ** 2

        elif total_range:
            # Returns the entire channel wavelength ranges effective area

            # iterate through lists for multiplication (otherwise type error occurs)
            reflectance = [i * j for i, j in zip(var['primary'], var['secondary'])]
            transmission_efficiency = [i * j for i, j in zip(var['fp_filter'], var['ent_filter'])]

            # equation:
            eff_area = (var['geoarea'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'contam'] * var['ccd']

            self.eff_area = eff_area
        else:
            # Returns effective area at the position of the center wavelength for the channel specified

            # # This test confirms that I need to index the values differently than just the index of the channel
            # for value in var['fp_filter']:
            #     if 0.348 < values < 0.349: # trying to replicate table data
            #         index = var['fp_filter'].index(value) # 684   != 94
            #         print(index, value)
            #         print(self.get_wavelength_range(channel)[index]) # 93.9 !!

            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    # print(n,  value)
                    index = n
                    break

            reflectance = var['primary'][index] * var['secondary'][index]
            transmission_efficiency = var['fp_filter'][index] * var['ent_filter'][index]

            # assert statement, these two should be the same
            # print(len(var['fp_filter']))
            # print(var['wavenumsteps'])

            # print what's being called, perform test if not close to what's expected
            # print(var['primary'][index], var['secondary'][index], var['fp_filter'][index], var['ent_filter'][index])
            # print(reflectance, transmission_efficiency)

            # equation:
            eff_area = (var['geoarea'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'contam'][index] * var['ccd'][index]
            self.eff_area = eff_area

            if print_comparison_ratio:
                compare_eff_area = [0.312, 1.172, 2.881, 1.188, 1.206, 0.063, 0.045, 0.0192, 0.0389][
                    self.channel_list.index(channel)]
                print('eff_area: ', compare_eff_area, eff_area, 'ratio of comparison: ', eff_area / compare_eff_area)

        return self.eff_area

    def wavelength_response(self, channel, compare_genx=False, compare_table=False, use_eff_area=False):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature.

        # formula from paper
        [R(wavelength, t) = A eff (wavelength, t)G(wavelength)]


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
                 count rate for spectrum?  : units cm^2 DN phot^-1
        """

        var = self.get_channel_data(channel)

        # units
        angstrom = u.angstrom
        # angstrom = 10**(-8) * u.cm
        dn_per_photon = u.count / u.photon



        # gain values from Boerner paper
        gain_table = {94: 2.128 * dn_per_photon, 131: 1.523 * dn_per_photon, 171: 1.168 * dn_per_photon,
                      195: 1.024 * dn_per_photon, 211: 0.946 * dn_per_photon,
                      304: 0.658 * dn_per_photon, 335: 0.596 * dn_per_photon, 1600: 0.125 * dn_per_photon,
                      1700: 0.118 * dn_per_photon}
        ccd_gain = {94: 18.3, 131: 17.6, 171: 17.7, 193: 18.3, 211: 18.3, 304: 18.3, 335: 17.6, 1600: 0.125,
                    1700: 0.118}
        read_noise_table = {94: 1.14, 131: 1.18,  171: 1.15, 193: 1.20, 211: 1.20, 304: 1.14, 335: 1.18, 1600: 1.15,
                    1700: 1.15}


        # gives instrument response values listed in Table 2
        if compare_genx:
            index = self.channel_list.index(channel)
            self.inst_response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][index]

        # Uses Gain from table Table 12 to calculate instrument response
        elif compare_table:
            for key, gain in gain_table.iteritems():
                if str(key) == str(channel):
                    inst_response = self.effective_area(key) * gain
                    # print(key, gain) # ~ 2.1
                    self.inst_response = inst_response

        elif use_eff_area:
            # TODO: TRY  gain = elecperphot / elecperdn --- no elecperphot in genx
            gain = var['elecperdn']# 18.3 for all channels from .genx file
            # print(gain)# ~ 18.3
            inst_response = self.effective_area(channel) * (gain  * (u.count / u.photon) )
            self.inst_response = inst_response


        # calculate the values
        else:
            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    # print('        ', n,  value)
                    index = n
                    break
                    # G(wavelength) = (12398 / wavelength / 3.65)g

            # the photon energy in eV
            ev = (12398.4953 * u.eV * u.angstrom * u.nm)
            wavelength =  (var['wave'][index] * u.nm )
            ev_per_wavelength = ev / wavelength
            print('elecperphot: ', ev_per_wavelength)
            #^^^ assert for wavelength 94, expect about 131.89 nm
            # these two should NOT be the same,
            # print('wave', var['wave'][index] , var['wave'][channel])


            # these are approximately the same
            calc_electron_per_ev =  (1* u.electron ) / (3.98 *  u.eV)
            electron_per_ev = var['elecperev'] * (u.electron / u.eV)
            print('e', electron_per_ev, calc_electron_per_ev)

            # rename this, more like electron per wavelength
            electron_to_photon = (ev_per_wavelength *  electron_per_ev)
            print('etp', electron_to_photon)

            electron_per_dn = var['elecperdn'] * (u.electron / u.count)
            platescale = var['platescale'] # very small...  10**(-11)

            # Gain(wavelength)  calulation
            calculated_gain = electron_to_photon  * (18.3 * (u.electron / u.count) * dn_per_photon)
            #  calculated_gain =  ( (12398.4953 / (float(channel) / 0.1)) * (1 /  3.98 ) )* (18.3) # channel* value *(10**(-5))


            # print('channel:', channel )
            # print('index:', index )
            # print('value:', value )
            # print('platescale: ', var['platescale'])
            # print('gain: ', [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][
            #     self.channel_list.index(channel)], calculated_gain * electron_per_dn)  # ~ 17 or 18.3 dn/phot

            eff_area = self.effective_area(channel, total_range=False, print_comparison_ratio=False)

            instr_response = eff_area * (calculated_gain / electron_per_dn)
            self.inst_response = instr_response

        return self.inst_response  # want units of cm**2 * count / photon

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
