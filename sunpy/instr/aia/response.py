from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"  # github: tdwilkinson

"""
This program performs characterization of the Atmospheric Imaging Assembly (AIA) instrument used on the Solar Dynamics
Observatory satellite. The wavelength and temperature response allow for photometric calibration using instrument
parameters. The instrument does not provide spectroscopic information, so there is no way to directly apply a
wavelength-dependent calibration to the data. These AIA response functions were calculated using AIA instrument
information obtained using Solar Soft Ware .genx files, and also integrating ChiantiPy ion temperature and density
information.

Utilizing this code one can:
    - Analyze AIA instrument response that tells efficiency of each channel using instrument properties.
    - Plot effective area vs wavelengths to see where on the instrument the response is peaked in wavelength.
    - Obtain temperature response based on chianti spectral contribution functions as a function of temperature.


"""

import numpy as np
import pandas as pd
import os

import astropy.units as u
import astropy.constants as constants

from scipy import integrate

from .aia_read_genx2table import aia_instr_properties_to_table
from .make_ion_contribution_table import save_contribution_csv


class Response():


    def __init__(self, channel_list, **kwargs):  # data_table,

        path_to_genx_dir = kwargs.get('path_to_genx_dir', '')
        version = kwargs.get('version', 6)

        self.data_table = aia_instr_properties_to_table(path_to_genx_dir, channel_list,
                                                        version, save=True)

        self.channel_list = channel_list

    def get_channel_data(self, channel):
        """
        This definition returns a specific row in self.data_table containing all of the channel data.
        Individual properties of the channel can be obtained by indexing the property.
        ie:  min_wavelength = get_channel_data(94)['wavemin']

        Parameters
        ---------
        :param channel: int
            specify which channel of which to obtain data


        Returns
        ------
        :return: Table row that acts like dictionary
        """

        # define the location of the channel data to return
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

            # returns  the table row that is a dictionary
            return channel_data

    def wavelength_range_index(self, channel):
        """
        Often times in ChiantiPy, the index of the wavelength is more useful than the the wavelength itself.
        This function returns the index of a channels center wavelength from the wavelength array.

        Parameter
        ---------
        :param channel: int
            specify which channel of which to obtain data

        Returns
        -------
        :return: the index value
        """

        var = self.get_channel_data(channel)
        num = var['number_wavelength_intervals']
        intervals = var['wavelength_interval']
        min_wave = var['minimum_wavelength']
        wave_range = np.arange(0, float(num)) * float(intervals) + float(min_wave)

        wavelength_range = var['wavelength_range']

        # assert wavelength_range.all() == wave_range.all() # check these are the same

        # obtain index for values as specific wavelength in wavelength range
        for n, value in enumerate(wavelength_range):
            if channel < value < channel + 1:
                index = n
                break

        return index

    def calculate_effective_area(self, channel,  use_genx_effarea=True):
        """
        AIA photometric calibration was obtained by making component-level measurements of all the optical elements in
        the AIA telescopes (mirrors, filters, and CCD), and combining those measurements analytically to produce a model
        of the system performance. This function contains information about the efficiency of the telescope optics.

        Formula:
        Effective Area = Geo_Area * Reflectance of Primary and Secondary Mirrors * Transmission Efficiency of filters
                        * Quantum Efficiency of CCD * Correction for time-varying contamination

        Parameters
        ----------
        :param channel, int
            the wavelength center of the channel to probe

        :param use_genx_effarea, book
            Default is to use the imported genx effective area, otherwise, the effective area will be caluculated using
            imported variables.

        Variables:
        ----------
        These are the variables needed to calculate the effective area for this instrument:

        Reflectance = The combined reflectance of the primary and secondary mirrors.
        Transmission_efficiency = The focal plane and entrance filters efficiency.
        Quantum Efficiency = The QE of the CCD
        Contamination = The correction for the time-varying contamination. This variable will change the most between
                        .genx versions

        Returns
        -------
        effective_area: np.array
            returns the effective area of the instrument over the entire wavelength range

        """

        var = self.get_channel_data(channel)

        if use_genx_effarea:
            # returns the channel wavelength effective area strait from the .genx files
            self.effective_area = var['effective_area'] * u.cm ** 2

        else:
            # Returns the entire channel wavelengths  effective area as calculated with .genx values

            # iterate through lists for multiplication (otherwise type error occurs)
            reflectance = [i * j for i, j in
                           zip(var['primary_mirror_reflectance'], var['secondary_mirror_reflectance'])]
            transmission_efficiency = [i * j for i, j in
                                       zip(var['focal_plane_filter_efficiency'], var['entire_filter_efficiency'])]

            # effective area equation:
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'] * var['quantum_efficiency_ccd']

            self.effective_area = eff_area




    def calculate_system_gain(self, channel, use_genx_values= True):
        """
        The CCD camera system gain is calculated using  a standard conversion of photons to detected electrons with the
        camera gain.

        Parameters
        --------
        :param: channel, int
            the wavelength center of the channel to probe

        :param: use_genx_values, bool
            Default is to use the values imported by genx files. Otherwise, the ccd system gain is calculated with the
            equation from Boerner et al 2012.

        """

        var = self.get_channel_data(channel)

        # conversions
        hc = (constants.h * constants.c).to(u.eV * u.angstrom)
        electron_per_ev = (1 * u.electron) / (3.98 * u.eV)

        if use_genx_values:
            # calculate the index of values to use from wavelength range
            channel_index = self.wavelength_range_index(channel)

            # elecperphot in genx starting with the photon energy in eV
            wavelength = (var['wavelength_range'][channel_index] * u.angstrom)
            ev_per_wavelength = hc / wavelength

            # converts energy of wavelength from eV to electrons
            genx_electron_per_ev = var['electron_per_ev'] * (u.count / u.eV)
            electron_per_wavelength = (ev_per_wavelength * genx_electron_per_ev)

            # gain = elecperphot / elecperdn
            self.system_gain = electron_per_wavelength / var['electron_per_dn']


        else:
            # calculate entire response using calculated eff_area and ccd gain from table12
            ccd_gain = {94: 18.3, 131: 17.6, 171: 17.7, 193: 18.3, 211: 18.3, 304: 18.3, 335: 17.6, 1600: 0.125,
                        1700: 0.118}

            # the photon energy in eV
            wavelength = (float(channel) * u.angstrom)
            ev_per_wavelength = hc / wavelength

            # converts energy of wavelength from eV to electrons
            electron_per_wavelength = (ev_per_wavelength * electron_per_ev)

            # gain of ccd from genx file:
            ccdgain = ccd_gain[channel] * (u.electron / u.count)

            # Gain(wavelength)  calculation
            calculated_gain = electron_per_wavelength * (u.count / u.photon)

            self.system_gain = calculated_gain / ccdgain




    def calculate_wavelength_response(self, channel, use_genx_values=True):
        """
        Describes the (AIA) instrument wavelength response by calculating effective area as a function of wavelength for
         the strongest emission lines present in the solar feature. This should display a peaked value around the
         channel wavelength centers.

        formula:
        Wavelength Response = Effective Area * Gain of the Intrument System

        Parameters
          --------
        :param: channel, int
            the wavelength center of the channel to probe

        :param: use_genx_values, bool
            Default is to use the values imported by genx files. Otherwise, the ccd system gain is calculated with the
            equation from Boerner et al 2012.



        Returns
        -------
        :return: response dictionary
            response dictionary contains the response per wavelength of effective area where the keys are the channel
            wavelength centers and both the responses and wavelenghts are returned as values.


        NOTE:
        The energy of a photon E = hf = hc / lambda . Expressing hc in units of eV = 12398.4953

        """

        # conversions
        dn_per_photon = u.count / u.photon

        var = self.get_channel_data(channel)

        # default is to use use effective area generated with a genx file
        if use_genx_values:

            self.calculate_system_gain(channel, use_genx_values=True)
            wavelength_response = var['effective_area'] * (self.system_gain * dn_per_photon)

        # otherwise an effective area calculation is completed using instrument properties from the genx file.
        else:
            self.calculate_system_gain(channel)
            self.calculate_effective_area(channel)
            wavelength_response = self.effective_area * self.system_gain

        wavelengths = np.arange(0, var['number_wavelength_intervals']) * var['wavelength_interval'] + var[
            'minimum_wavelength']

        out_response = {}
        out_response[channel] = {'response': wavelength_response, 'wavelength': wavelengths}

        self.wavelength_response = out_response




    def get_wavelength_response_functions(self):
        """
        Calc wavelength response functions for each channel in self.channel_list

        :return: dictionary
            each key is the channel and the value is the wavelength response array
        """

        out_dictionary = {}

        for channel_wavelength in self.channel_list:
            self.calculate_wavelength_response(channel_wavelength, use_genx_values=True)
            out_dictionary[channel_wavelength] = self.wavelength_response[channel_wavelength]

        self.wavelength_response_functions = out_dictionary




    def calculate_temperature_response(self, channel, trapz1=False, trapz2=False, **kwargs):
        """
        Calculates temperature for various features using ChiantiPy which analyze the sensitivity of light from the
        plasma per temperature as measured by AIA. The instrument temperature-response functions also provide some basic
        insight into the interpretation of the images from the various channels.

        Parameters
        ----------
        :param: channel, int
            The wavelength center of the channel to probe

        :param: trapz1, bool
            An alternative integration method for calculating this response. The default is simpson integration.

        :param: trapz2, bool
            An alternative integration method for calculating this response. the default is simpson integration.

        :param: temperature_array, array
            The default to this 10. ** (5.0 + 0.05 * np.arange(61.)) to give an outgoing array that goes from 10**5 to
            10**8 in temperature range. This can be adjusted to be smaller and become less time intensive for the ion
            calculation.

        :return: dict
            temperature response - dictionary with channel wavelength center as keys and temperatures and temperature
            response values as items

        """
        print(channel)

        # isobaric model pressure is 10^15 - needed?
        pressure = 10 ** 15 * (u.Kelvin / u.cm ** 3)

        # default temperature
        temperature_array = kwargs.get('temperature', 10. ** (5.0 + 0.05 * np.arange(61.)))

        # GENERALIZATION: try density = t * pressure to get array ( / temperature * kb)
        density = 1.e+9

        # calculate wavelength response # expect this to be run first
        self.calculate_wavelength_response(channel)

        # save ion contribution into a dataframe for faster computation
        ion_contribution_dataframe = kwargs.get('ion_contribution_dataframe', 'test.csv')

        # check for path. does this work with defining your own ion data frame ^^ ?
        if not os.path.exists(ion_contribution_dataframe):
            # load wavelength response function for all channels
            self.get_wavelength_response_functions()

            print('WARNING: Making an Ion datatable will take a while!')
            save_contribution_csv(self.channel_list, self.wavelength_response_functions, temperature_array,
                                  density=density)

        data = pd.read_csv('test.csv', delimiter=',', header=0, index_col=0)

        # from saved data table, define all wavelengths in range around channel into 2D matrix
        array = []
        discrete_wavelengths = []

        # from saved data table, only search wavelength columns around the channel in question
        for number in np.arange(channel - 3, channel + 3, .001):
            if str(number) in data.keys():
                store_per_wavelength = []
                discrete_wavelengths.append(number)
                for ion_contribution_per_temp in data[str(number)]:
                    store_per_wavelength.append(ion_contribution_per_temp)
                array.append(store_per_wavelength)

        # change the vertical values to be horizontal before interpolation
        matrix = np.vstack((array)).T

        # evaluate wavelength response function at specific wavelengths, DOES affect magnitude:
        response_interpolated = []
        for wavelength in discrete_wavelengths:
            response_interpolated.append(self.wavelength_response[channel]['wavelength'][wavelength])

        # multiply ion contributions array by wave-length response array
        response = matrix * response_interpolated

        # integrate through wavelength range to sample specifically near discrete_wavelength values
        # method 1: composite trapezoidal integration
        if trapz1:
            temp_response = integrate.cumtrapz(response, discrete_wavelengths)

        # method 2: another trapezoidal integration
        elif trapz2:
            temp_response = integrate.trapz(response, discrete_wavelengths)
        else:
            # method 3: simpson integration
            temp_response = integrate.simps(response, discrete_wavelengths)

        # define temperature response dictionary
        self.temperature_response = {'temperature': temperature_array, 'temperature response': temp_response}




    def get_temperature_response_functions(self):
        '''

        Calculate the  temperature response functions for each channel in self.channel_list

        :return: dictionary
            Each key is the channel and the value is a dictionary with temperatures and temperature response array
        '''



        out_dictionary = {}

        for channel_wavelength in self.channel_list:
            self.calculate_temperature_response(channel_wavelength)
            out_dictionary[channel_wavelength] = self.temperature_response[channel_wavelength]

        self.temperature_response_functions = out_dictionary
