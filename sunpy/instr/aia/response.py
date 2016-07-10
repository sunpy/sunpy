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

    def __init__(self, channel_list, **kwargs):  # data_table,

        path_to_genx_dir = kwargs.get('path_to_genx_dir', '')
        version = kwargs.get('version', 6)

        # If no data table exists, make the data table
        # if not os.path.exists('channel_properties_' + str(version) + '.csv'):
        #     aia_read_genx2table.aia_instr_properties_to_table(path_to_genx_dir, channel_list, properties, version, save=True)

        # load in astropy table  ( quicker than reading genx files each time)
        # self.data_table = Table.read('channel_properties_' + str(version) + '.csv')

        self.data_table = aia_read_genx2table.aia_instr_properties_to_table(path_to_genx_dir, channel_list,
                                                                            version, save=True)


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

        num = self.get_channel_data(channel)['number_wavelength_intervals']
        intervals = self.get_channel_data(channel)['wavelength_interval']
        min_wave = self.get_channel_data(channel)['minimum_wavelength']

        if channel == 'all':
            # wavelength_range = np.arange(0, float(self.data_dictionary[0][numwavesteps])) * float(
            #     (self.data_dictionary[0][wave_intervals])) + float(self.data_dictionary[0][min_wavelength])
            pass
        else:
            wavelength_range = np.arange(0, float(num)) * float(intervals) + float(min_wave)
            # assert should match var['wavelength_range']

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
            self.eff_area = var['effective_area'] * u.cm ** 2

        elif total_range:
            # Returns the entire channel wavelength ranges effective area

            # iterate through lists for multiplication (otherwise type error occurs)
            reflectance = [i * j for i, j in zip(var['primary_mirror_reflectance'], var['secondary_mirror_reflectance'])]
            transmission_efficiency = [i * j for i, j in zip(var['focal_plane_filter_efficiency'], var['entire_filter_efficiency'])]

            # equation:
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'] * var['quantum_efficiency_ccd']

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

            reflectance = var['primary_mirror_reflectance'][index] * var['secondary_mirror_reflectance'][index]
            transmission_efficiency = var['focal_plane_filter_efficiency'][index] * var['entire_filter_efficiency'][index]

            # assert statement, these two should be the same
            # print(len(var['fp_filter']))
            # print(var['wavenumsteps'])

            # print what's being called, perform test if not close to what's expected
            # print(var['primary'][index], var['secondary'][index], var['fp_filter'][index], var['ent_filter'][index])
            # print(reflectance, transmission_efficiency)

            # equation:
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'][index] * var['quantum_efficiency_ccd'][index]
            self.eff_area = eff_area

            if print_comparison_ratio:
                compare_eff_area = [0.312, 1.172, 2.881, 1.188, 1.206, 0.063, 0.045, 0.0192, 0.0389][
                    self.channel_list.index(channel)]
                print('eff_area: ', compare_eff_area, eff_area, 'ratio of comparison: ', eff_area / compare_eff_area)

        return self.eff_area

    def get_wavelength_response(self, channel, use_response_table2=False, use_calc_effarea_table2_gain=False,
                                use_genx_values=False):
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
        dn_per_photon = u.count / u.photon

        # gain values from Boerner paper
        gain_table2 = {94: 2.128, 131: 1.523, 171: 1.168, 193: 1.024, 211: 0.946, 304: 0.658, 335: 0.596, 1600: 0.125,
                       1700: 0.118}
        ccd_gain = {94: 18.3, 131: 17.6, 171: 17.7, 193: 18.3, 211: 18.3, 304: 18.3, 335: 17.6, 1600: 0.125,
                    1700: 0.118}

        # returns the instrument response values listed in Table 2
        if use_response_table2:
            index = self.channel_list.index(channel)
            response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][index]
            self.wavelength_response = response
        # Uses Gain from table Table 12 to calculate instrument response
        elif use_calc_effarea_table2_gain:
            gain = gain_table2[channel]
            # for key, gain in gain_table2.iteritems():
            #     if str(key) == str(channel):
            # print(key, gain) # ~ 2.1
            self.wavelength_response = self.effective_area(channel) * gain * dn_per_photon

        elif use_genx_values:
            # no elecperphot in genx - have to make it with imported properties

            # calculate the index of values to use from wavelength range
            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    # print('        ', n,  value)
                    index = n
                    break

            # the photon energy in eV
            ev = (12398.4953 * u.eV * u.angstrom)
            wavelength = (var['wavelength_range'][index] * u.angstrom)
            ev_per_wavelength = ev / wavelength

            electron_per_ev = var['electron_per_ev'] * (u.count / u.eV)

            # converts energy of wavelength from eV to electrons
            electron_per_wavelength = (ev_per_wavelength * electron_per_ev)

            # gain = elecperphot / elecperdn
            gain = electron_per_wavelength / var['electron_per_dn']
            # print(gain)# 18.3 for all channels from .genx file

            self.wavelength_response = var['effective_area'][index] * (gain * dn_per_photon)


        else:
            # calculate entire response using calculated eff_area and ccd gain from table12

            # the photon energy in eV
            ev = (12398.4953 * u.eV * u.angstrom)
            wavelength = (float(channel) * u.angstrom)
            ev_per_wavelength = ev / wavelength
            # print('elecperphot: ', ev_per_wavelength)
            # ^^^ assert for wavelength 94, expect about 131.89 nm
            # these two should NOT be the same,
            # print('wave', var['wave'][index] , var['wave'][channel])

            calc_electron_per_ev = (1 * u.electron) / (3.98 * u.eV)
            # print('e', electron_per_ev, calc_electron_per_ev)            # these are approximately the same

            # converts energy of wavelength from eV to electrons
            electron_per_wavelength = (ev_per_wavelength * calc_electron_per_ev)
            # print('epl', electron_per_wavelength)

            # gain of ccd from genx file:
            ccdgain = ccd_gain[channel] * (u.electron / u.count)
            # print(ccdgain)

            # Gain(wavelength)  calulation
            calculated_gain = electron_per_wavelength * dn_per_photon

            # print('channel:', channel )
            # print('index:', index )
            # print('value:', value )
            # print('platescale: ', var['platescale'])
            # print('gain: ', [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][
            #     self.channel_list.index(channel)], calculated_gain * electron_per_dn)  # ~ 17 or 18.3 dn/phot

            calc_eff_area = self.effective_area(channel, total_range=False, print_comparison_ratio=False)

            self.wavelength_response = calc_eff_area * (calculated_gain / ccdgain)

        return self.wavelength_response  # want units of cm**2 * count / photon

    def emissivity():
        """

        - Use  chianti.core.mspectrum.intensityratio to get spectral line intensities and show relative emissivity per line.
            chianti.core.continuum has places to access freebound and freefree emissions
        - create a line list  to show line intensities

        input:


        :return: generates chianti object - continuum - chianti model with line list
        """

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
