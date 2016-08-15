from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""
The Atmospheric Imaging Assembly (AIA) instrument used on the Solar Dynamics Observatory satellite do not provide spectroscopic information, so there is no way to directly apply a wavelength-dependent calibration to the data.
These AIA response functions were calculated using AIA instruement information obtained using Solar Soft Ware .genx files, and also integrating ChiantiPy ion temperature and density information.

Utilizing this code one can:
    - Analyze AIA instrument response that tells efficiency of each channel using instrument properties.
    - Plot effective area vs wavelengths to see where on the instrument the response is peaked in wavelength.
    - Obtain temperature response based on chianti spectral contribution functions as a function of temperature/density.


- Test with sunpy/instr/tests/test_aia.py

An example of how to run this code:

"insert example text here!"

"""

import numpy as np
import pandas as pd
import os

import astropy.units as u
import astropy.constants as constants

from scipy import integrate
from scipy.interpolate import interp1d

from .aia_read_genx2table import aia_instr_properties_to_table
from .make_ion_contribution_table import save_contribution_csv
import sunpy.data.test as test


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

        self.data_table = aia_instr_properties_to_table(path_to_genx_dir, channel_list,
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

    def wavelength_range_index(self, channel):
        """
        wavelength range is calculated here for plotting becomes more general after calling a specific channel
        - needed because r.get_channel_data(94)['wavelength_range'] is an array full of the number 25.0
        :param channel:
        :return:
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

    def calculate_effective_area(self, channel, total_range=False, use_genx_effarea=False,
                                 print_comparison_ratio=False):
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

        if use_genx_effarea:
            self.effective_area = var['effective_area'] * u.cm ** 2

        elif total_range:
            # Returns the entire channel wavelength ranges effective area

            # iterate through lists for multiplication (otherwise type error occurs)
            reflectance = [i * j for i, j in
                           zip(var['primary_mirror_reflectance'], var['secondary_mirror_reflectance'])]
            transmission_efficiency = [i * j for i, j in
                                       zip(var['focal_plane_filter_efficiency'], var['entire_filter_efficiency'])]

            # equation:
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'] * var['quantum_efficiency_ccd']

            self.effective_area = eff_area

        else:
            # Returns effective area float value at the position of the center wavelength for the channel specified
            # want to move away from this and just return an array with zeros everywhere but where the channel is specified
            channel_index = self.wavelength_range_index(channel)

            # import instrument properties
            reflectance = var['primary_mirror_reflectance'][channel_index] * var['secondary_mirror_reflectance'][
                channel_index]
            transmission_efficiency = var['focal_plane_filter_efficiency'][channel_index] * \
                                      var['entire_filter_efficiency'][
                                          channel_index]

            # equation calculation :
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'][channel_index] * var['quantum_efficiency_ccd'][channel_index]
            self.effective_area = eff_area

            if print_comparison_ratio:
                # used for testing - move to test_aia?
                compare_eff_area = [0.312, 1.172, 2.881, 1.188, 1.206, 0.063, 0.045, 0.0192, 0.0389][
                    self.channel_list.index(channel)]
                print('eff_area: ', compare_eff_area, eff_area, 'ratio of comparison: ', eff_area / compare_eff_area)

    def calculate_system_gain(self, channel, use_table2_gain=False, use_genx_values=False):
        """

        :return:
        """

        var = self.get_channel_data(channel)

        # conversions
        hc = (constants.h * constants.c).to(u.eV * u.angstrom)
        electron_per_ev = (1 * u.electron) / (3.98 * u.eV)

        if use_table2_gain:
            # Uses Gain from table Table 12 to calculate instrument response
            # gain values from Boerner paper
            gain_table2 = {94: 2.128, 131: 1.523, 171: 1.168, 193: 1.024, 211: 0.946, 304: 0.658, 335: 0.596,
                           1600: 0.125, 1700: 0.118}
            self.system_gain = gain_table2[channel]


        elif use_genx_values:
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

    def calculate_wavelength_response(self, channel, use_genx_values= True,  **kwargs):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature.

        # formula from paper
        [R(wavelength, t) = A eff (wavelength, t)G(wavelength)]


        Parameters
        ----------
        :keyword


        Returns
        -------
        :return: float, array describing the response per wavelength of effective area (wavelength response)
                 count rate for spectrum?  : units cm^2 DN phot^-1

        NOTE:

        The energy of a photon E = hf = hc / lambda . Expressing hc in units of eV = 12398.4953

        """

        # conversions
        dn_per_photon = u.count / u.photon

        var = self.get_channel_data(channel)

        # # move to test
        # if use_table2_gain:
        #     self.calculate_effective_area(channel)
        #     self.calculate_system_gain(channel, use_table2_gain=True)
        #     self.wavelength_response = self.effective_area * self.system_gain * dn_per_photon

        # # move to test
        # elif use_response_table2:
        #     # returns the instrument response values listed in Table 2
        #     i = self.channel_list.index(channel)
        #     response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][i]
        #     self.wavelength_response = response

        # default is to use use effective area generated with a genx file
        if use_genx_values:

            self.calculate_system_gain(channel, use_genx_values=True)
            wavelength_response = var['effective_area'] * (self.system_gain * dn_per_photon)

        # otherwise an effective area calculation is completed using instrument properties from the genx file.
        else:
            self.calculate_system_gain(channel)
            self.calculate_effective_area(channel, total_range=True)
            wavelength_response = self.effective_area * self.system_gain

        wavelengths = np.arange(0, var['number_wavelength_intervals'])*var['wavelength_interval']+var['minimum_wavelength']

        out_response = {}
        out_response[channel] = {'response': wavelength_response, 'wavelength': wavelengths}

        self.wavelength_response = out_response




    def get_wavelength_response_functions(self):
        """
        calc wavelength response functions for each channel in self.channel_list

        :return: dictionary
            each key is the channel and the value is the wavelength response array
        """
        out_dictionary = {}

        # if array
        for channel_wavelength in self.channel_list:
            out_dictionary[channel_wavelength] = self.calculate_wavelength_response(channel_wavelength, use_genx_values = True)

        self.wavelength_response_functions = out_dictionary






    def calculate_temperature_response(self, channel, trapz1=False, trapz2=False, **kwargs):
        """
        calculates temperature for various features using ChiantiPy

        instrument temperature-response functions provides some basic insight into the interpretation of
                the images from the various channels.

        input: string, data file
        channel: string, indicates channels used

        :return: dict

            dictionary with temperatures and temperature response as items

        """

        # isobaric model pressure is 10^15 - needed?
        pressure = 10 ** 15 * (u.Kelvin / u.cm ** 3)

        # default temperature
        temperature_array = kwargs.get('temperature', 10. ** (5.0 + 0.05 * np.arange(61.)))

        # TODO: try density = t * pressure to get array ( / temperature * kb)
        density = 1.e+9
        # print('t: ', np.log(t))
        # print('t',t)

        # calculate wavelength response # expect this to be run first
        self.calculate_wavelength_response(channel)


        # save ion contribution into a dataframe for faster computation
        ion_contribution_dataframe = kwargs.get('ion_contribution_dataframe', 'test_ion_contributions.csv')

        # check for path. does this work with defining your own ion data frame ^^ ?
        if not os.path.exists(ion_contribution_dataframe):
            print('WARNING: Making an Ion datatable will take a while!')
            save_contribution_csv(self.channel_list, temperature_array, density = density)

        data = pd.read_csv('test_ion_contributions.csv', delimiter = ',', header = 0, index_col = 0)
        # data = pd.read_csv('top_ions_contributions.csv', delimiter = ',', header = 0, index_col = 0)
        # data = pd.read_csv('all_ions_contributions.csv', delimiter = ',', header = 0, index_col = 0)




        # from saved dataframe, put values near channel into 2D matrix
        array = []
        discrete_wavelengths = []

        for number in np.arange(channel - 3, channel + 3, .001): # auto-sort?
            if str(number) in data.keys():
                store_per_wavelength = []
                discrete_wavelengths.append(number)
                for ion_contribution_per_temp in data[str(number)]:
                    store_per_wavelength.append(ion_contribution_per_temp)
                array.append(store_per_wavelength)
        # print('array', array, len(array))
        matrix = np.vstack((array)).T

        # print('matrix:' , len(matrix)) # 61
        # print(discrete_wavelengths) # returns array   (1,3) for channel 94, only three ions probed - check

        # check that these are the same length -- good! 8751 == 8751
        # print( len(self.wavelength_response[channel]['wavelength']), len(self.wavelength_response[channel]['response']))


        # Wavelength Response Interpolation
        # here more ions is good
        # response_interpolated =  interp1d(discrete_wavelengths,  self.wavelength_response[channel]['wavelength'], self.wavelength_response[channel]['response'])[0]
        # response_interpolated = np.interp(discrete_wavelengths,  self.wavelength_response[channel]['wavelength'], self.wavelength_response[channel]['response'])[0]

        # evaluate wavelength response function at specific wavelengths, DOES affect magnitude: (-10**2 from above interpolation)
        response_interpolated = []
        for wavelength in discrete_wavelengths:
            response_interpolated.append(self.wavelength_response[channel]['wavelength'][wavelength])

        # # multiply ion contributions array by wave-length response array
        response = matrix * response_interpolated


        # move to test: check: these are the same length because wavelength range in defining response functions
        # print('response', response, len(response), len(response[0]))
        # print(len(wrange))

        # integrate through wavelength range, default = simpsons, sample specifically near discrete_wavelength values
        # method 1: composite trapezoidal integration
        if trapz1:
            temp_response = integrate.cumtrapz(response, discrete_wavelengths)

        # method 2:
        elif trapz2:
            temp_response = integrate.trapz(response, discrete_wavelengths)
        else:
            # method 3: simpson integration
            temp_response = integrate.simps(response, discrete_wavelengths)


        # define temperature response dictionary
        # print(len(temperature_array), len(temp_response), print(temp_response))
        self.temperature_response = {'temperature': temperature_array, 'temperature response': temp_response}




    def get_temperature_response_functions(self):
        '''

        calc temperature response functions for each channel in self.channel_list

        :return: each key is the channel and the value is a dictionary with temperatures, density, and temperature response array
        '''

        out_dictionary = {}

        # if array
        for channel_wavelength in self.channel_list:
            out_dictionary[channel_wavelength] = self.calculate_temperature_response(channel_wavelength)

        self.temperature_response_functions = out_dictionary

