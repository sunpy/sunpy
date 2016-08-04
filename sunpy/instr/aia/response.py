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

import astropy.units as u
import astropy.constants as constants
import chianti.core as ch
import chianti.constants as c
from scipy import integrate

from .aia_read_genx2table import aia_instr_properties_to_table
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

    def calculate_wavelength_response(self, channel, band_width=2.5,
                                      use_response_table2=False, use_table2_gain=False, use_genx_values=True, **kwargs):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature.

        # formula from paper
        [R(wavelength, t) = A eff (wavelength, t)G(wavelength)]


        Parameters
        ----------
        :keyword

        band_width in angstroms - range around channel wavelength that is probed, all other values return as zeros in wavelength response array

        cropped = test case of returning a cropped array
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

        # default is to use these values
        if use_genx_values:

            self.calculate_system_gain(channel, use_genx_values=True)
            wavelength_response = var['effective_area'] * (self.system_gain * dn_per_photon)


            # if total range removed. float value returned, moving away from this!
            # else:
            #     self.wavelength_response = var['effective_area'][channel_index] * (
            #         self.system_gain * dn_per_photon)

        else:
            self.calculate_system_gain(channel)
            self.calculate_effective_area(channel, total_range=True)
            wavelength_response = self.effective_area * self.system_gain



            # wavelength range selected, conversion to index by multiplying by ten.
            # This is because the wavelength array is of dispersed wavelengths.
            # # + 100 in index = + 10 wavelengths in array
            # # + 10 in index = + 1 wavelength in array
        self.range = int(band_width * 10)
        channel_index = self.wavelength_range_index(channel)

        # mask the wavelength array with zeros except in a small wavelength range near channel
        for index, value in enumerate(wavelength_response):
            if index - self.range < channel_index < index + self.range:
                continue
            else:
                wavelength_response[index] = 0.0

        self.wavelength_response = wavelength_response

    def calculate_contribution_function(self, channel, ion_object, chianti_Gofnt=False,
                                        intensity=False):
        """
        information on plasma and atomic physics governing how matter emits at a given temperature

            - Chianti contribution function tells how much a transition contributes to the emissivity at a given
                temperature or line intensity per unit emission measurement
        input:

        :return:
        """
        # default temperature - make this global?
        # t = 10. ** (5.5 + 0.05 * np.arange(21.))
        # print('t: ', t)
        var = self.get_channel_data(channel)

        # wavelength range for chianti ion object needs to be given in form [start, stop]

        # method 1: using the whole wavelength range given from SSW
        # start_wvl = var['wavelength_range'][0]
        # end_wvl = var['wavelength_range'][-1]

        # method 2: using only a wavelength band around the channel ie. channel 94, using 92.5-95.5

        wrange = var['wavelength_range']
        channel_index = self.wavelength_range_index(channel)
        start_wvl = wrange[channel_index - self.range]  # default set in wavelength range to 2.5 angstroms
        end_wvl = wrange[channel_index + self.range]
        wavelength_range = [start_wvl, end_wvl]

        if chianti_Gofnt:
            #   refactored chianti gofnt
            # print(
            #     'ERROR: chianti.core.ion.gofnt does not return self.Gofnt until after prompts and plots for each ion.')

            # Chianti Contribution Function
            ion_object.gofnt(wvlRange=wavelength_range, top=1, verbose=0, plot=False)

            self.temperatures = ion_object.Gofnt['temperature']
            self.contribution_function = ion_object.Gofnt['gofnt']

            # print('len gofnt', len(ion_object.Gofnt['gofnt']))
            # print('eDensity', ion_object.Gofnt['eDensity'])
            # print('index', ion_object.Gofnt['index'])
            # print('wvl', ion_object.Gofnt['wvl'])

            # self.Gofnt = {'temperature': outTemperature, 'eDensity': outDensity, 'gofnt': gofnt, 'index': g_line,
            #               'wvl': wvl[g_line]}

    def get_ion_contribution_array(self, channel, temperature, top_ions=True, all_chianti_ions=False,
                                   solar_chromosphere_ions=False, test=False, **kwargs):
        """

        :return:
        """
        # looping through elements and ions:
        ion_array = []

        # test case: store iron ion names
        if test:
            # for level in range(5, 28, 1):
            #     ion_array.append('fe' + '_' + str(level))

            # just iron in 94 give a pretty good result for channel 94 only
            ion_array = ['fe_10', 'fe_18', 'fe_20']

        if top_ions:
            # default: top ions selected from chianti linelist that have high intensities in near each channel wavelength
            ion_array = ['al_10', 'ca_14', 'ca_15', 'ca_17', 'fe_3', 'fe_7', 'fe_9', 'fe_11', 'fe_12', 'fe_14', 'fe_16',
                         'fe_18', 'fe_19', 'fe_20', 'fe_21', 'fe_23', 'ni_14', 'ni_20', ]

        # # test case: store ion names of elements with high solar abundance - basically same result as running all ions
        # # (based on logarithmic abundance chart, higher than ~ 7 dex)
        if solar_chromosphere_ions:
            # solar_elements = ['h', 'he', 'c', 'n', 'o', 'ne', 'si', 's', 'ar', 'ca', 'ni', 'fe']
            # for element in solar_elements:
            #     for level in range(1, 18, 1): # max 38
            #         ion_array.append(element + '_' + str(level))
            # or
            #
            # store every ion from Felding and Widing (1993) abundances: bibcode 1995ApJ...443..416L
            ion_array = ['o_5', 'o_6', 'ne_8', 'mg_7', 'mg_10', 'si_6', 'si_7', 'si_8', 'si_9', 'si_10', 's_8', 's_9',
                         'fe_14', 'fe_15', 's_10', 's_11', 's_12', 's_13', 'fe_8', 'fe_9', 'fe_10', 'fe_11', 'fe_12',
                         'fe_13', 'ni_11', 'ni_12', 'ni_13', 'ni_17', 'ni_18']

        # # test case: store every element in chianti ion names - takes ~ 10 mins per channel, ~ 80 mins for all channels
        if all_chianti_ions:
            for element in c.El:
                for level in range(1, 38, 1):
                    ion_array.append(element + '_' + str(level))

        ion_contribution_dataframe = kwargs.get('ion_contribution_dataframe', 'contributions.csv')
        # if ion_contribution_dataframe:
        #     # read contributions from dataframe
        #     # move to another definition?
        #     # df = pd.dataframe()
        #     #ion_contributions = array
        #     pass
        # else:

        # find contribution function for all ions
        ion_contributions = []
        for i in ion_array:
            try:
                ion_object = ch.ion(i, temperature=temperature, eDensity=1.e+9, em=1.e+27)
                self.calculate_contribution_function(channel, ion_object, chianti_Gofnt=True)
                ion_contributions.append(self.contribution_function)
                print(i)
            except (AttributeError, KeyError, UnboundLocalError, IndexError):
                # doesn't save ions if no lines found in selected interval
                #                      missing information (Elvlc file missing)
                #                      cause errors in ChiantyPy (Zion2Filename, self.Ip)
                # print(i, ' contribution not calculated')
                pass

        # sum all contributions for all ions
        self.contributions_array = sum(ion_contributions)

        # move to test: check: confirms sum did the job!
        # print('ca: ', contributions_array)
        # print('ic: ', ion_contributions, len(ion_contributions))

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
        t = kwargs.get('temperature', 10. ** (5.0 + 0.05 * np.arange(61.)))
        # print('t: ', np.log(t))
        # print('t',t)

        # calculate wavelength response # expect this to be run first
        self.calculate_wavelength_response(channel, cropped=False)

        # get ion contributions
        self.get_ion_contribution_array(channel, temperature=t)

        # multiply ion contributions array by wave-length response array
        response = []
        for ion_contribution_per_temp in self.contributions_array:
            response.append(ion_contribution_per_temp * self.wavelength_response)

        # move to test: check: these are the same length because wavelength range in defining response functions
        # print('save', response, len(response), len(response[0]))
        # print(len(wrange))

        # integrate through wavelength range, default = simpsons
        # method 1: composite trapezoidal integration
        if trapz1:
            temp_response = integrate.cumtrapz(response, axis=1)

        # method 2:
        elif trapz2:
            temp_response = integrate.trapz(response, axis=1)
        else:
            # method 3: simpson integration
            temp_response = integrate.simps(response, axis=1)

        # test if multiplying by pressure and/or plate scale give better accuracy
        # temp_response = tr #* ( 1 / pressure) * (1 / platescale)
        # print(temp_response, len(temp_response))

        # define temperature response dictionary
        self.temperature_response = {'temperature': self.temperatures, 'temperature response': temp_response}
