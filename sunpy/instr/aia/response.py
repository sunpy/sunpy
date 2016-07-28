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
from scipy import integrate

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

    def get_wavelength_range_channel_index(self, channel):
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

        self.channel_index = index


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
            self.get_wavelength_range_channel_index(channel)

            # import instrument properties
            reflectance = var['primary_mirror_reflectance'][self.channel_index] * var['secondary_mirror_reflectance'][self.channel_index]
            transmission_efficiency = var['focal_plane_filter_efficiency'][self.channel_index] * var['entire_filter_efficiency'][
                self.channel_index]

            # equation calculation :
            eff_area = (var['geometric_area_ccd'] * u.cm ** 2) * reflectance * transmission_efficiency * var[
                'ccd_contamination'][self.channel_index] * var['quantum_efficiency_ccd'][self.channel_index]
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
            self.get_wavelength_range_channel_index(channel)

            # elecperphot in genx starting with the photon energy in eV
            wavelength = (var['wavelength_range'][self.channel_index] * u.angstrom)
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

    def calculate_wavelength_response(self, channel, total_range=False, use_response_table2=False,
                                      use_table2_gain=False, use_genx_values=False, **kwargs):
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
            i = self.channel_list.index(channel)
            response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][i]
            self.wavelength_response = response

        elif use_genx_values:

            self.get_wavelength_range_channel_index(channel)
            wavelength_range = kwargs.get('wavelength_range', [self.channel_index - 25, self.channel_index + 25])

            self.calculate_system_gain(channel, use_genx_values=True)
            if total_range:
                self.wavelength_response = var['effective_area'] * (self.system_gain * dn_per_photon)
            elif wavelength_range:
                self.wavelength_response = var['effective_area'][wavelength_range] * (self.system_gain * dn_per_photon)
            else:
                self.wavelength_response = var['effective_area'][self.channel_index] * (self.system_gain * dn_per_photon)

        else:
            self.calculate_system_gain(channel)
            self.calculate_effective_area(channel, total_range=True)

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

        self.emissivity_wavelength_array = fe14.Emiss['wvl']  # 41745 length
        self.emissivity = fe14.Emiss['emiss']  # 21 values of 41745 length

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

    def get_most_intense_line(self, channel, wavelength_range, ion_object):
        """

        :return:
        """
        temperatures = ion_object.Temperature
        emissivity, emissivity_wavelength_array = self.get_emissivity(channel, wavelength_range, ion_object)

        # TODO: want to use ionWeb.gofntSelectLines to get self.toplines and self.wvlchoices
        # couldn't figure out how to access it...

        # # most intense lines  direct from chianti GOFNT
            # return topline_index, top_lines
        wvl = ion_object.Intensity["wvl"]
        igvl = range(len(wvl))
        nlines = len(igvl)
        igvl = np.take(igvl, wvl[igvl].argsort())
        # find the top most intense lines
        # if top > nlines: # default 10
        #     top = nlines
        #
        intensity = ion_object.Intensity['intensity']
        maxIntens = np.zeros(nlines, 'Float64')
        for iline in range(nlines):
            maxIntens[iline] = intensity[:, igvl[iline]].max()
        for iline in range(nlines):
            if maxIntens[iline] == maxIntens.max():
                maxAll = intensity[:, igvl[iline]]
        igvlsort = np.take(igvl, np.argsort(maxIntens))
        # print('igvlsort = ', igvlsort)
        topLines = igvlsort[-1:]  # previously top, default 10 in code
        # print(' topLines = ', topLines)
        maxWvl = '%5.3f' % wvl[topLines[-1]]
        #        maxline=topLines[-1]
        #
        topLines = topLines[wvl[topLines].argsort()]

        return topLines, maxWvl

    def calculate_contribution_function(self, channel, wavelength_range, ion_object, chianti_Gofnt=False,
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
        var = self.get_channel_data(channel)

        if chianti_Gofnt:
            #   refactored chianti gofnt
            # print(
            #     'ERROR: chianti.core.ion.gofnt does not return self.Gofnt until after prompts and plots for each ion.')

            # Chianti Contribution Function
            ion_object.gofnt(wvlRange=wavelength_range, top=1, verbose=0, plot=False)
            gofnt = ion_object.Gofnt['gofnt']

            self.temperatures = ion_object.Gofnt['temperature']
            self.contribution_function = gofnt

            # self.Gofnt = {'temperature': outTemperature, 'eDensity': outDensity, 'gofnt': gofnt, 'index': g_line,
            #               'wvl': wvl[g_line]}

        elif intensity:

            # intensity INCLUDES elemental abundance and ionization fraction
            #           can it be used to calculate gofnt without plotting
            #           as recommended by kdere of chiantipy
            #           fe14.Intensity.keys() = ['em', 'intensity', 'wvl', 'lvl2', 'lvl1',
            #                        'pretty2', 'ionS', 'avalue', 'pretty1', 'integrated', 'obs']

            ion_object.intensity(wvlRange=wavelength_range)

            emissivity_array = ion_object.Intensity['em']  # array of one value em=1.e+27 given in def of ion
            # ^^^^ EMISSIVITY IS USED IN SSW ^^^^

            # # Attempt to get top line
            # peak_line, maxwvl = self.get_most_intense_line(channel, wavelength_range, ion_object)
            # print('top line: ', peak_line, maxwvl)



            # find index in wavelength range for the line at peak value
            # index = np.argmin(np.abs(ion_object.Intensity['wvl'] - float(maxwvl)))
            # print('iw:', intensity_wavelengths[index], index)
            # print(intensity_array[:,index])
            index = ion_object.Intensity['intensity'].index(ion_object.Intensity['intensityTop'])
            print('index: ', index)

            # Intensity['intensity'] is double nested
            gofnt = ion_object.Intensity['intensity'][:, index] / ion_object.EDensity
            print('gofnt: ', gofnt)

            self.temperatures = ion_object.Temperature
            self.contribution_function = gofnt


        else:

            # not sure we need a model spectrum? just in case...
            # spec = self.spectrum(temperature=t, wavelength_range=wavelength_range, specific_ion=specific_ion, channel=channel)

            # emissivity defined separate from Contribution Function (Gofnt)
            self.get_emissivity(channel=channel, wavelength_range=wavelength_range, specific_ion=ion_object)

            # ssw multiplies by platescale
            platescale = self.get_channel_data(channel)['plate_scale']
            print('platescale', platescale)  # check that out...

            # top intensity lines in range defined separate from contribution function
            ion_object.intensity(wvlRange=wavelength_range)
            top_line_index, nlines = self.get_most_intense_line(channel=channel, wavelength_range=wavelength_range,
                                                                ion_object=ion_object)

            # Contribution function variables defined alternately
            ion_object.populate(temperature=t)  # creates fe14.Population containing the level population information
            edensity = ion_object.Population['eDensity']
            ion_eq = ion_object.IoneqOne  # the ionization equilibrium for the selected ion as a function of temperature.
            g_ioneq = ion_eq / edensity  # divided by the electron density
            abundance = ion_object.Abundance

            # gofnt function per line
            self.contribution_function = abundance * g_ioneq * self.emissivity[top_line_index]
            # print('gofnt: ', self.contribution_function)




    def make_ion_contribution_array(self, channel, temperature, wavelength_range, test = False, all_chianti_ions = False, solar_chromosphere_ions = False ):
        """

        :return:
        """
        # looping through elements and ions:
        ion_array = []

        # test case: store iron ion names
        if test:
            # ion_array =  [ 'fe_14', 'fe_15']
            for level in range(1, 15, 1):
                ion_array.append('fe' + '_' + str(level))

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

            ion_array = ['o_5', 'o_6', 'ne_8', 'mg_7', 'mg_10', 'si_6', 'si_7', 'si_8', 'si_9', 'si_10', 's_8', 's_9', 'fe_14',
                          'fe_15', 's_10', 's_11', 's_12', 's_13', 'fe_8', 'fe_9', 'fe_10', 'fe_11', 'fe_12', 'fe_13','ni_11', 'ni_12', 'ni_13', 'ni_17', 'ni_18']


        # # test case: store every element in chianti ion names - takes ~ 10 mins per channel, ~ 80 mins for all channels
        if all_chianti_ions:
            for element in c.El:
                for level in range(1, 38, 1):
                    ion_array.append(element + '_' + str(level))
            # # TODO: create data structure with contribution structure of all elements for easy access


        # find contribution function for all ions
        ion_contributions = []
        for i in ion_array:
            try:
                ion_object = ch.ion(i, temperature=temperature, eDensity=1.e+9, em=1.e+27)
                self.calculate_contribution_function(channel, wavelength_range, ion_object, chianti_Gofnt=True)
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

        # check: confirms sum did the job!
        # print('ca: ', contributions_array)
        # print('ic: ', ion_contributions, len(ion_contributions))

    def calculate_temperature_response(self, channel, simp = False, trapz1 = False, trapz2 = False, **kwargs):
        """
        calculates temperature for various features using ChiantiPy

        instrument temperature-response functions provides some basic insight into the interpretation of
                the images from the various channels.

        input: string, data file
        channel: string, indicates channels used

        :return: dict

            dictionary with temperatures and temperature response as items

        """
        var = self.get_channel_data(channel)




        # temperature response variables: wavelengthrange needs to be given as [start, stop]
        wrange = var['wavelength_range']
        # print('wave', wrange)

        start_wvl = self.get_channel_data(channel)['minimum_wavelength']
        end_wvl = wrange[-1]
        wavelength_range = [start_wvl, end_wvl]
        print('wavelength range1: ', wavelength_range)

        #
        start_wvl = channel - 30
        end_wvl = channel + 30
        wavelength_range = [start_wvl, end_wvl]
        print('wavelength range2: ', wavelength_range)

        #
        # obtain index for values as specific wavelength in wavelength range
        self.get_wavelength_range_channel_index(channel)

        start_wvl = wrange[self.channel_index - 25]
        end_wvl = wrange[self.channel_index + 25]
        wavelength_range = [start_wvl, end_wvl]
        print('wavelength range3: ', wavelength_range)





        # isobaric model pressure is 10^15 - needed?
        pressure = 10 ** 15 * (u.Kelvin / u.cm ** 3)

        # default temprature
        t = kwargs.get('temperature', 10. ** (5.5 + 0.05 * np.arange(21.)))
        # print('t',t)

        # # test with wavelength response as full array?
        self.calculate_wavelength_response(channel, total_range = True, use_genx_values=True)

        platescale = var['plate_scale']  # * (1 / u.second)

        # get ion contributions
        self.make_ion_contribution_array(channel, temperature = t, wavelength_range = wavelength_range, all_chianti_ions=True)


        # multiply ion contributions array by wave-length response array
        response = []
        for ion_contribution_per_temp in self.contributions_array:
            response.append(ion_contribution_per_temp * self.wavelength_response)

        # check: these are the same length because wavelength range in defining response functions
        print('save', response, len(response), len(response[0]))
        print(len(wrange))

        # integrate through wavelength range
        # method 1: simpson integration
        if simp:
            temp_response = integrate.simps(response, axis = 1)

        # method 2: composite trapezoidal integration
        if trapz1:
            temp_response = integrate.cumtrapz(response, axis = 1)

        # method 3:
        if trapz2:
            # transverse = response.transpose()
            # resp = np.vstack((response))
            # print('r', resp, len(resp), len(resp[0]) #  2D 21 (,8761
            temp_response = integrate.trapz(response, axis = 1)

        # method 4: general purpose integration
        # # define function for integrate     # 2D axis, integrate along wavelength axis
        # def temperature_response_function(temp):
        #     return # not sure
        #
        # tr, err = integrate.quad(temperature_response_function(response, wavelength_range), start_wvl, end_wvl)
        # temp_response.append(tr)

        # test if multiplying by pressure and/or plate scale give better accuracy
        # temp_response = tr #* ( 1 / pressure) * (1 / platescale)
        # print(temp_response, len(temp_response))

        # define temperature response dictionary
        self.temperature_response = {'temperature': self.temperatures, 'temperature response': temp_response}
