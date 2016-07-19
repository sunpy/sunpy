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

import sunpy
import sunpy.data.test as test
from .aia_read_genx2table import aia_instr_properties_to_table

# running github ChiantiPy!
import chianti.core as ch
import chianti.constants as c
import matplotlib.pyplot as plt


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

    def get_wavelength_range(self, channel):
        """
        wavelength range is calculated here for plotting becomes more general after calling a specific channel
        - needed because r.get_channel_data(94)['wavelength_range'] is an array full of the number 25.0
        :param channel:
        :return:
        """

        num = self.get_channel_data(channel)['number_wavelength_intervals']
        intervals = self.get_channel_data(channel)['wavelength_interval']
        min_wave = self.get_channel_data(channel)['minimum_wavelength']

        if channel == 'all':
            # TODO:  return an array / ndarray of wavelength data
            # wavelength_range = np.arange(0, float(self.data_dictionary[0][numwavesteps])) * float(
            #     (self.data_dictionary[0][wave_intervals])) + float(self.data_dictionary[0][min_wavelength])
            pass
        else:
            wavelength_range = np.arange(0, float(num)) * float(intervals) + float(min_wave)

            return wavelength_range

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
            # Returns effective area at the position of the center wavelength for the channel specified

            # find the index in arrays from ssw for desired wavelength
            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    # print(n,  value)
                    index = n
                    break

            # import instrument properties
            reflectance = var['primary_mirror_reflectance'][index] * var['secondary_mirror_reflectance'][index]
            transmission_efficiency = var['focal_plane_filter_efficiency'][index] * var['entire_filter_efficiency'][
                index]

            # equation calculation :
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'][index] * var['quantum_efficiency_ccd'][index]
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

            # obtain index for values as specific wavelength
            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    index = n
                    break

            # elecperphot in genx starting with the photon energy in eV
            wavelength = (var['wavelength_range'][index] * u.angstrom)
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
            calculated_gain = electron_per_wavelength * dn_per_photon

            self.system_gain = calculated_gain / ccdgain

    def calculate_wavelength_response(self, channel, use_response_table2=False, use_table2_gain=False,
                                      use_genx_values=False):
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

        if use_table2_gain:
            self.calculate_effective_area(channel)
            self.calculate_system_gain(channel, use_table2_gain=True)
            self.wavelength_response = self.effective_area * self.system_gain * dn_per_photon


        elif use_response_table2:
            # returns the instrument response values listed in Table 2
            index = self.channel_list.index(channel)
            response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][index]
            self.wavelength_response = response

        elif use_genx_values:
            # obtain index for values as specific wavelength in wavelength range
            for n, value in enumerate(self.get_wavelength_range(channel)):
                if channel < value < channel + 1:
                    index = n
                    break

            self.calculate_system_gain(channel, use_genx_values=True)
            self.wavelength_response = var['effective_area'][index] * (self.system_gain * dn_per_photon)

        else:
            self.calculate_system_gain(channel)
            self.calculate_effective_area(channel)

            # TODO: fix error in units count**2?!   # want units of cm**2 * count / photon
            self.wavelength_response = self.effective_area * self.system_gain

    def get_emissivity(self, channel, wavelength_range, specific_ion):
        """
        Find the emissivity for each line
         - develop an emissivity (chianti) spectral structure based on plasma  properties (line intensities?/ emissivity?)

        :return: generates chianti object - continuum - chianti model with line list

        Notes:


        - create a line list  to show line intensities


        # solar abundances in chianti!
        chianti abundfile  default = sun_photospheric_1998_grevesse which includes the abundances of
                                    Grevesse and Sauval, 1998, Space Science Reviews, 85, 161.
                                    Boerner paper uses: Feldman and Widing (1993) ?


        # ionization equilibrium in chianti!
        chianti ioneqfile = The default value is chianti which includes the ionization equilibrium calculations
                            of Dere, et al., 2009, Astronomy and Astrophysics, 498, 915 and are considered to be based
                            on the best ionization and recombination rates currently available.
        """

        t = 10. ** (5.5 + 0.05 * np.arange(21.))
        fe14 = ch.ion('fe_14', temperature=t, eDensity=1.e+9, em=1.e+27)

        # ion needs form [start, stop] -- what's the best method for choosing the wavelength range?
        wrange = self.get_wavelength_range(channel)
        start_wvl = self.get_channel_data(channel)['minimum_wavelength']
        end_wvl = wrange[-1]
        wavelength_range = [start_wvl, end_wvl]

        # # define for emissitivity function attributes - not required ?
        # fe14.populate(temperature=t)  # creates fe14.Population containing the level population information
        # fe14.Population['population'] # level populations for specific ion
        # fe14.IoneqOne # the ionization equilibrium for the selected ion as a function of temperature.

        # create fe14.Emiss containing the line emissivity information
        fe14.emiss(wvlRange=wavelength_range)

        emissivity_wavelength_array = fe14.Emiss['wvl']  # 41745 length
        emissivity = fe14.Emiss['emiss']  # 21 values of 41745 length

        return emissivity, emissivity_wavelength_array

    def spectrum(self, temperature, wavelength_range, fe_14, channel):
        """
        generate a spectral model for various solar features to be used in modules
                 idl uses this to make a 'structure' object...
                 maybe I can use it to make a line list?

        - want to read values of physical parameters to determine feature
            chianti.core.continuum has places to access freebound and freefree emissions

        input:

        :return:
        """

        # # continuum to get freebound and freefree emissions
        # # want:  ioneq and abundance to default values
        # fe14_cont = ch.continuum('fe_14', temperature=temperature)
        # fe14_cont.freeFreeEmiss(wvl=wavelength_range)
        # fe14_cont.freeBoundEmiss(wvl=channel)
        #
        # fe14_cont.freeFree(wvl=channel)
        # fe14_cont.freeBound(wvl=channel)
        # # fe14_cont.IoneqName  = chianti ?
        #
        #
        # # use spectrum to get synthetic spectrum with minwl, maxwl, pressure, goft, ioneqfile  of ions?
        # # ie) ch.spectrum(temperature, eDensity, wavelength, filter, label, elementList, ionList, minAbund, doContinuum, em, )
        # # ie ) ch.spectrum(temperature=t, eDensity=fe14.Population['eDensity'], wavelength=channel )
        #
        # # test:  attribute error: continuum instance has no attribute 'FreeFree' in ch.spectrum()
        # # fe14.spectrum() # creates fe14.Spectrum contain the line and continuum spectrum information
        # fe14.spectrum(wvl=wavelength_range, em=emissivity_array, eDensity=fe14.Population['eDensity'])
        #
        # return spectrum

    def get_most_intense_line(self, channel, wavelength_range, specific_ion):
        """

        :return:
        """
        temperatures = specific_ion.Temperature
        emissivity, emissivity_wavelength_array = self.get_emissivity(channel, wavelength_range, specific_ion)

        # TODO: want to use ionWeb.gofntSelectLines to get self.toplines and self.wvlchoices
        # couldn't figure out how to access it...

        # most intense lines  direct from chianti GOFNT
        igvl = range(len(emissivity_wavelength_array))
        nlines = len(emissivity_wavelength_array)
        maxEmiss = np.zeros(nlines, 'Float64')
        for iline in range(nlines):
            maxEmiss[iline] = emissivity[igvl[iline]].max()
        for iline in range(nlines):
            if maxEmiss[iline] >= maxEmiss.max():
                maxAll = emissivity[igvl[iline]]  # used in defining output
        igvlsort = np.take(igvl, np.argsort(maxEmiss))
        topLines = igvlsort[-nlines:]

        # print(topLines)
        top_lines = topLines[emissivity_wavelength_array[topLines].argsort()]  # better way to do this?

        # should plot just like GOFNT, test:
        # out = []
        # for iline in range(top):
        #     tline = topLines[iline]
        #     out.append(emissivity[tline]/maxAll)
        #
        #     plt.loglog(temperatures, emissivity[tline]/maxAll)
        # plt.show()

        # index position of highest intensity line
        topline_index = emissivity_wavelength_array[top_lines[0]]

        return topline_index, nlines

    def calculate_contribution_function(self, channel, wavelength_range, specific_ion, chianti_Gofnt=False,
                                        intensity=False):
        """
        information on plasma and atomic physics governing how matter emits at a given temperature

            - Chianti contribution function tells how much a transition contributes to the emissivity at a given
                temperature or line intensity per unit emission measurement
        input:

        :return:
        """
        # default temperature - make this global?
        t = 10. ** (5.5 + 0.05 * np.arange(21.))

        if chianti_Gofnt:
            #   refactor chianti gofnt? time?
            print(
                'ERROR: chianti.core.ion.gofnt does not return self.Gofnt until after prompts and plots for each ion.')

            # Chianti Contribution Function  ## PROMPTS AND PLOT ###
            specific_ion.gofnt(wvlRange=wavelength_range, top=1, verbose=0)
            gofnt = specific_ion.Gofnt['gofnt']
            self.contribution_function = gofnt

            # self.Gofnt = {'temperature': outTemperature, 'eDensity': outDensity, 'gofnt': gofnt, 'index': g_line,
            #               'wvl': wvl[g_line]}

        elif intensity:

            # --- [WIP] --- if possible/worth it?
            # intensity may be used for plotting one line (alt. to gofnt)
            # AND it INCLUDES elemental abundance and ionization fraction.
            # can it be used to calculate these values instead of using gofnt?

            specific_ion.intensity(wvlRange=wavelength_range)

            # fe14.Intensity contain the line intensities information
            # size 21 from len temperatures - each array then is len 41745 - same as gofnt output
            intensity_array = specific_ion.Intensity['intensity']  # array of [[values]]!
            emissivity_array = specific_ion.Intensity['em']  # array of one value em=1.e+27 given in def of ion
            print('ratio:', intensity_array / emissivity_array)

            # intensity array or line intensity?
            self.contribution_function = intensity_array / emissivity_array


        else:

            # not sure we need a model spectrum? just in case...
            # spec = self.spectrum(temperature=t, wavelength_range=wavelength_range, specific_ion=specific_ion, channel=channel)

            # emissivity defined separate from Contribution Function (Gofnt)
            emissivity, em_wavelength_array = self.get_emissivity(channel=channel, wavelength_range=wavelength_range,
                                                                  specific_ion=specific_ion)

            # ssw multiplies by platescale
            platescale = self.get_channel_data(channel)['plate_scale']

            # top intensity lines in range defined separate from contribution function
            top_line_index, nlines = self.get_most_intense_line(channel=channel, wavelength_range=wavelength_range,
                                                                specific_ion=specific_ion)

            # Contribution function variables defined alternately
            specific_ion.populate(temperature=t)  # creates fe14.Population containing the level population information
            edensity = specific_ion.Population['eDensity']
            ion_eq = specific_ion.IoneqOne  # the ionization equilibrium for the selected ion as a function of temperature.
            g_ioneq = ion_eq / edensity  # divided by the electron density
            abundance = specific_ion.Abundance

            # gofnt function per line
            self.contribution_function = abundance * g_ioneq * emissivity[top_line_index]
            print('gofnt: ', self.contribution_function)

    def calculate_temperature_response(self, channel):
        """
        calculates temperature for various features using ChiantiPy

        instrument temperature-response functions provides some basic insight into the interpretation of
                the images from the various channels.

        input: string, data file
        channel: string, indicates channels used

        :return: dict

            dictionary with temperatures and temperature response as items

        """

        # load in wavelength response
        self.calculate_wavelength_response(channel, use_genx_values=True)

        # all elements we can check... need all of them?
        ion_array = []
        for element in c.El:
            for level in range(1, 38, 1):
                ion_array.append(element + '_' + str(level))

        # temperature response variables:

        # wavelengthrange needs to be given as [start, stop] -    probably a better way to et this...
        wrange = self.get_wavelength_range(channel)
        start_wvl = self.get_channel_data(channel)['minimum_wavelength']
        end_wvl = wrange[-1]
        wavelength_range = [start_wvl, end_wvl]
        # print('wavelength range: ', wavelength_range)

        # isobaric model pressure is 10^15
        pressure = 10 ** 15 * (u.Kelvin / u.cm ** 3)

        # default temperature
        t = 10. ** (5.5 + 0.05 * np.arange(21.))

        # define ion
        fe14 = ch.ion('fe_14', temperature=t, eDensity=1.e+9, em=1.e+27)
        # define ion temperatures
        temperatures = fe14.Temperature  # same as fe14.Population['temperature']  # length 21

        # pulled this to separate function until debugged. remerge as needed.
        self.calculate_contribution_function(channel, wavelength_range, fe14)

        temp_response = self.contribution_function * self.wavelength_response
        print('temp response: ', temp_response)

        # TODO: loop through array of ions (from chianti elements), Too many, need to choose some. But still Testing:
        # gofnt_array = []
        # for i in ion_array:
        #     print(i)   # breaks on h_2 .. :( 'ion' has no attribute 'Nlvls' for h_2
        #     specific_ion = ch.ion( i , temperature=t, eDensity=1.e+9, em=1.e+27)
        #     temperatures = specific_ion.Temperature  # same as fe14.Population['temperature']  # length 21
        #     contribution_function = self.get_contribution_function(channel, wavelength_range, specific_ion)
        #     gofnt_array.append(contribution_function)
        # temp_response = gofnt_array * wavelength_response
        # print('temp response: ', self.temperature_response)

        self.temperature_response =  {'temperature': temperatures, 'temperature response': temp_response}
