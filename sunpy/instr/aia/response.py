"""
Provides functions for calculating wavelength and temperature response of
SDO/AIA
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.constants as constants
from scipy import integrate

from .aia_read_genx2table import aia_instr_properties_to_table
from .make_ion_contribution_table import save_contribution_csv

__author__ = ["Tessa D. Wilkinson","Will Barnes"]


class Response():
    """
    This class calculates the wavelength and temperature response functions for the 7 EUV channels of the Atmospheric Imaging Assembly (AIA) instrument onboard the Solar Dynamics Observatory spacecraft. These AIA response functions were calculated using AIA instrument information obtained from  Solar Software (SSW) .genx files. Temperature response functions are calculated using spectroscopic information from the CHIANTI atomic database and the ChiantiPy package.

    Parameters
    ----------

    Examples
    --------

    Notes
    -----

    References
    ----------
    """


    channel_colors = {94:'#ff3a3a',131:'#6060ff',171:'#f1de1f',193:'#4cec4c',211:'#ed64c6',
                        335:'#45deed',304:'k',1600:'b',1700:'g',4500:'r'}

    def __init__(self, channel_list=[94,131,171,193,335,211,304], ssw_path='', version=6):
        tmp = os.path.join(ssw_path,'sdo','aia','response','aia_V{}_{}_fullinst.genx')
        instrument_files = [tmp.format(version,'all'),tmp.format(version,'fuv')]
        self._get_channel_info(channel_list,instrument_files)

    def _get_channel_info(self,channel_list,instrument_files):
        """
        Get instrument info for channels in channel list. Creates self._channel_info

        Notes
        -----
        This will probably change once instrument data is acquired through other means.
        """
        data_table = aia_instr_properties_to_table(channel_list,instrument_files)
        self._channel_info = {}
        for c in channel_list:
            index = channel_list.index(c)
            self._channel_info[c] = {property : data_table[property][index] for property in data_table.columns.keys()}

    def calculate_effective_area(self, channel):
        """
        AIA photometric calibration was obtained by making component-level measurements of all the optical elements in
        the AIA telescopes (mirrors, filters, and CCD), and combining those measurements analytically to produce a model
        of the system performance. This function contains information about the efficiency of the telescope optics.

        Formula:
        Effective Area = Geo_Area * Reflectance of Primary and Secondary Mirrors * Transmission Efficiency of filters
                        * Quantum Efficiency of CCD * Correction for time-varying contamination

        These are the variables needed to calculate the effective area for this instrument:

        Reflectance = The combined reflectance of the primary and secondary mirrors.
        Transmission_efficiency = The focal plane and entrance filters efficiency.
        Quantum Efficiency = The QE of the CCD
        Contamination = The correction for the time-varying contamination. This variable will change the most between
                        .genx versions

        Parameters
        ----------
        channel : int
            the wavelength center of the channel to probe

        Returns
        -------
        effective_area: np.array
            returns the effective area of the instrument over the entire wavelength range

        """
        reflectance = self._channel_info[channel]['primary_mirror_reflectance']*self._channel_info[channel]['secondary_mirror_reflectance']
        transmission_efficiency = self._channel_info[channel]['focal_plane_filter_efficiency']*self._channel_info[channel]['entrance_filter_efficiency']

        # effective area equation:
        effective_area = self._channel_info[channel]['geometric_area_ccd']*reflectance*transmission_efficiency*self._channel_info[channel]['ccd_contamination']*self._channel_info[channel]['quantum_efficiency_ccd']

        return effective_area

    def calculate_system_gain(self, channel):
        """
        The CCD camera system gain is calculated using  a standard conversion of photons to detected electrons with the camera gain.

        Parameters
        --------
        channel : int
            the wavelength center of the channel

        Notes
        -----
        The energy of a photon E = hf = hc / lambda . 12398.4953 = hc in eV
        Angstrom
        """
        # convert hc to convenient units
        hc = (constants.h * constants.c).to(u.eV * u.angstrom)
        #energy per photon in units of eV
        ev_per_photon = hc/(self._channel_info[channel]['wavelength'])/u.photon
        electron_per_ev = self._channel_info[channel]['electron_per_ev']
        electron_per_photon = ev_per_photon*electron_per_ev
        # gain = elecperphot / elecperdn
        return electron_per_photon/self._channel_info[channel]['electron_per_dn']

    def calculate_wavelength_response(self):
        """
        Describes the (AIA) instrument wavelength response by calculating effective area as a function of wavelength for the strongest emission lines present in the solar feature. This should display a peaked value around the channel wavelength centers.

        formula:
        Wavelength Response = Effective Area * Gain of the Intrument System

        Notes
        -----
        Does the platescale need to be included in this calculation?
        """

        self.wavelength_response = {}
        for channel in self._channel_info:
            system_gain = self.calculate_system_gain(channel)
            effective_area = self.calculate_effective_area(channel)
            self.wavelength_response[channel] = {'wavelength':self._channel_info[channel]['wavelength'], 'response':system_gain*effective_area}

    def peek_wavelength_response(self):
        """
        Quick plot of wavelength response functions versus wavelength
        """
        if not hasattr(self,'wavelength_response'):
            self.calculate_wavelength_response()

        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()
        for key in self.wavelength_response:
            ax.plot(self.wavelength_response[key]['wavelength'], self.wavelength_response[key]['response']/np.max(self.wavelength_response[key]['response']), label='{0} {1:latex}'.format(key,u.angstrom), color=self.channel_colors[key])

        ax.set_xlabel(r'$\lambda$ ({0:latex})'.format(self.wavelength_response[key]['wavelength'].unit))
        ax.set_ylabel(r'$R_i(\lambda)$ (norm.) ({0:latex})'.format(self.wavelength_response[key]['response'].unit))
        lambda_min = np.min(self.wavelength_response.keys())
        lambda_max = np.max(self.wavelength_response.keys())
        lambda_diff= np.fabs(lambda_max-lambda_min)*0.25
        ax.set_xlim([lambda_min-lambda_diff,lambda_max+lambda_diff])
        ax.legend(loc='best')
        plt.show()

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

        num = self._channel_info[channel]['number_wavelength_intervals']
        intervals = self._channel_info[channel]['wavelength_interval']
        min_wave = self._channel_info[channel]['minimum_wavelength']
        wave_range = np.arange(0, float(num)) * float(intervals) + float(min_wave)

        wavelength_range = self._channel_info[channel]['wavelength_range']

        # assert wavelength_range.all() == wave_range.all() # check these are the same

        # obtain index for values as specific wavelength in wavelength range
        for n, value in enumerate(wavelength_range):
            if channel < value < channel + 1:
                index = n
                break

        return index

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
