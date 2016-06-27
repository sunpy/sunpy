from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""
AIA response functions by integrating ChiantiPy

    - analyze AIA instrument response that tells efficiency of channel using instrument properties
    - then output a wavelength response
    - use spectral model (synthetic spectra) - to get ion emissivities as a function of temp/density
    - obtain temperature response based on chianti spectral contribution functions (emissivity)

AIA ccd's do not provide spectroscopic information, so there is no way to directly apply a wavlength-dependent calibration to the data.

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
import os

import sunpy
import sunpy.data.test as test


class Response():
    """
    class description here!

    Parameters
    ----------
    path_to_genx_file: string
        give a path to a instrument file from .genx output to initialize a dataframe for the class
    save_to_outfile: bool
        loads channel properties into a dataframe and saves to a .csv file
    csv: string
        alternately, give a path to the csv file to read in data information where the keys are channel wavelength centers and the values are properties

    Attributes
    ----------
    properties: array of strings
        The names of properties that are imported from SSW
    wavelength_centers: array of integers
        The center wavelength of each wavelength range or channel imaged with AIA
    dataframe: pandas dataframe
        The property information stored after being imported from SSW
    wavelength_range: array
        the range of wavelengths imaged per channel with AIA


    # Example:
    # >>> from sunpy.instr.aia import Response
    #
    # # load in dataframe to be used
    # >>> r = Response(path_to_genx_file='/media/gsoc/ssw_aia_response_data/aia_V6_all_fullinst', save_outfile=True)
    #
    # # call functions
    # >>> r.get_properties_per_channel(94)


    Notes:
    Currently the instr directory contains the file aia.py (previously aiaprep.py - I think?)
    This branch creates instr/aia directory and will have both this program response.py and aiaprep.py

    Feedback/ Thoughts are always welcome! -Thanks!
    """

    def __init__(self, path_to_genx_dir, wavelength, **kwargs):

        # wavelength center list and parameter list are used when defining dataframe
        self.wavelength_centers = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]

        self.properties = ['wave', 'effarea', 'geoarea', 'platescale',
                           'numfilters', 'wavemin', 'wavestep', 'wavenumsteps', 'wavelog',
                           'usecontam', 'contamthick', 'elecperdn',
                           'useerror', 'fp_filter', 'ent_filter', 'primary', 'secondary', 'ccd', 'contam',
                           'cross_area']  # 'usephottoelec', 'usephottodn',  'elecperdn','elecperev', not in 1600, 1700

        # keyword arguements
        overwrite = kwargs.get('overwrite', False)
        version = kwargs.get('version', 6)
        csv = kwargs.get('csv', 'channel_properties_' + str(version) + '.csv')

        # Load in the filename and define dataframe
        if os.path.exists(csv):
            if overwrite:
                os.remove(csv)
                print('overwriting: ' + 'channel_properties_' + str(version) + '.csv')
                self.aia_inst_genx_to_dataframe(path_to_genx_dir, version=version)
            else:
                print('using: ' + 'channel_properties_' + str(version) + '.csv')
                self.dataframe = pd.DataFrame.from_csv('channel_properties_' + str(version) + '.csv', header=0,
                                                       index_col=0, sep=',')

        else:
            print('making: ', 'channel_properties_' + str(version) + '.csv')
            self.aia_inst_genx_to_dataframe(path_to_genx_dir, version=version)

        # defines properties for one channel with wavelength keyword. implement all wavelengths at once?
        # def __call__(self, **kwargs): # possible way to implement wavelength after initial^^^
        # currently to plot, have to loop through array of wavelengths

        assert wavelength in self.wavelength_centers, "Choose one integer value from wavelength centers: ' + str(self.wavelength_centers) or 'all'"
        self.wavelength = wavelength

        print(self.dataframe[str(wavelength)]['ccd']) # could replace property_per_channel_dict with this implemented here
        var = self.property_per_channel_dict(wavelength)


        self.platescale = kwargs.get('platescale', var['platescale'])
        self.elecperdn = kwargs.get('elecperdn', var['elecperdn'])
        self.wavelength_range = kwargs.get('wave', var['wave'])
        self.wavemin = kwargs.get('wavemin', var['wavemin'])
        self.wavenumsteps = kwargs.get('wavenumsteps', var['wavenumsteps'])
        self.wavestep = kwargs.get('wavestep', var['wavestep'])
        # test: wavelength_range == wavemin+wavestep * wavenumsteps
        self.geometric_area = kwargs.get('geoarea', var['geoarea'])
        self.focal_plane_filter = kwargs.get('fp_filter', var['fp_filter'])
        self.entire_filter = kwargs.get('ent_filter', var['ent_filter'])
        self.primary = kwargs.get('primary', var['primary'])
        self.secondary = kwargs.get('secondary', var['secondary'])
        self.ccd = kwargs.get('ccd', var['ccd'])
        self.instrument_contamination = kwargs.get('contam', var['contam'])  # nan in 1600 and 1700
        # self.cross_area = kwargs.get('cross_area', var['cross_area'][0])
        # 1600, 1700 channels have arrays for some attributes #TODO: figure out why/ which is best
        if type(self.primary) == np.ndarray: # previously self.primary[0]
            self.primary = self.primary[0]
            self.secondary = self.secondary[0]

        # eff_area is not in all channels
        try:
            self.eff_area = kwargs.get('eff_area', var['effarea'])
        except KeyError:
            print(wavelength, 'has no eff_area keyword.')
            pass

    def aia_inst_genx_to_dataframe(self, input_directory, version):
        """

            Give the directory to a .genx directory to recursive search inside and return the the EUV and UV instrument files.
            output = list of paths to instrument files

        Parameters
        ----------
        path_to_file : string, the path location that leads to ssw_aia_response_data/aia_V6_all_fullinst.

        :param save: bool, if True, will save a dataframe containing all of the extracted instrument data to

        """
        check = []  # make sure not duplicated
        file_paths = []

        for root, dirs, files in os.walk(input_directory, topdown=False):

            for name in files:
                if (os.path.join(root, name)) not in check and name.endswith('_fullinst') and str(version) in str(name):
                    check.append(os.path.join(root, name))
                    file_paths.append(os.path.join(root, name))

        assert len(file_paths) != 0, 'Did not find aia_' + str(version) + '_fullinst files in directory.'

        # fill dataframe with desired properties and save csv file
        self.dataframe = self.save_properties_to_dataframe(file_paths, version)

    def save_properties_to_dataframe(self, files, version):
        """
        This definition reads the instrument file aia_V6_all_fullinst, which was obtained from ssw_aia_response_data inside
        ssw_aia_response_genx.tar.gz (an output file saved from SolarSoft Ware (SSW)). It will extract the instrument data and save
        it a dataframe file for easier access.


        Returns
        -------
        :returns self.dataframe: dataframe, stores the dictionary for easy access. Each column lists the channel information and each row is a different peroperty from aia_inst_genx.

        :returns 'channel_properties.csv', dataframe with wavelength centers for each channel is a column and listed in the column are the properties from the is returned where it has keys of the properties of the channel

        Notes:
        np.recarray store information with shape(1,0) and are quite nested.

        """

        # store in dataframe
        df = pd.DataFrame()

        for instr_file in files:
            # access np.recarray from .genx file
            ssw_array = readsav(instr_file)
            data = ssw_array['data']

            # pick out instrument files inside np.recarray
            for name in data.dtype.names:

                # Read in data from Channels from  A##_FULL  which match format in UV and EUV files
                if name.startswith('A') and name.endswith('_FULL') and name.find('THICK') < 0:

                    # target number in filename that matches wavelength_center
                    start = name.find('A')
                    end = name.find('_F')
                    wavelength = name[start + 1:end]

                    # UV files have 'thick' files that need to be filtered out for no duplicates
                    if int(wavelength) in self.wavelength_centers:

                        # match the properties in the rec.array to desired properties
                        variables = set(data[name][0].dtype.fields.keys()) & set(self.properties)

                        # then iterate through the properties to fill a dataframe with rec.array information
                        for property in variables:
                            try:
                                df.loc[property, wavelength] = data[name][0][property][0]
                            # contam and crossarea are sequences -- # TODO: work on this try statement
                            except ValueError:
                                df.loc[property, wavelength] = data[name][0][property]
                    else:
                        pass

        assert len(df) != 0, 'Data Frame is not loading from file.'

        #  save to .csv outfile
        df.to_csv('channel_properties_' + str(version) + '.csv')

        return df

    def property_per_channel_dict(self, wavelength_center=94):
        """
        Returns a dictionary for the channel whose wavelength center is inputed.

        Parameters
        ----------
        wavelength_center: int
            has to be one of the wavelengths in self.wavelength_centers


        Returns
        -------
        :returns A dictionary with the following keys and values:

        wave: array
            wavelengths imaged by the instrument and centered around the channel wavelength
        geoarea: float
            the geometric area of the each channel: EUV = 83.0 cm**2, UV = 30.8 cm**2
        ent_filter: array
            transmission efficiency of the entrance filter
        fp_filter: array of floats
            transmission efficiency of the focal plane filter
        contam: array of floats
            correction to adjust for the time varying contamination of the components
        secondary: array of floats
            reflectance of secondary mirror
        primary: array
            reflectance of the primary mirror
        ccd: array, floats
            the quantum efficiency of the ccd
        """
        assert wavelength_center in self.wavelength_centers, 'Choose one integer value from wavelength centers: ' + str(
            self.wavelength_centers)

        dict = {}

        # at target wavelength - pull out respective properties in dataframe and apply units
        for wavelength, values in self.dataframe.iteritems():
            if int(wavelength) == int(wavelength_center):
                dict['wavemin'] = values['wavemin'] * u.angstrom
                dict['effarea'] = values['effarea'] * u.cm ** 2
                dict['geoarea'] = values['geoarea'] * u.cm ** 2
                dict['platescale'] = values['platescale']
                dict['numfilters'] = values['numfilters']
                dict['wavestep'] = values['wavestep']
                dict['wavenumsteps'] = values['wavenumsteps']
                dict['wavelog'] = values['wavelog']
                dict['usecontam'] = values['usecontam']
                dict['contamthick'] = values['contamthick']
                dict['fp_filter'] = values['fp_filter']
                dict['ent_filter'] = values['ent_filter']
                dict['primary'] = values['primary']
                dict['secondary'] = values['secondary']
                dict['ccd'] = values['ccd']
                dict['contam'] = values['contam']
                dict['cross_area'] = values['cross_area']
                dict['wave'] = values['wave'] * u.angstrom

                # not in 1600 and 1700 channels
                try:
                    dict['elecperdn'] = values['elecperdn']  # loaded for now for comparison
                    dict['elecperev'] = values['elecperev']  # loaded for now for comparison
                except KeyError:
                    pass

                assert len(dict) != 0

                return dict

    def effective_area(self, compare_genx=False, compare_table=False):
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
        if compare_genx:
            self.eff_area = self.eff_area
        elif compare_table:
            self.eff_area = [0.312, 1.172, 2.881, 1.188, 1.206, 0.063, 0.045, 0.0192, 0.0389][
                self.wavelength_centers.index(self.wavelength)]
        else:
            # variables:
            wavelength = self.wavelength
            reflectance = self.primary * self.secondary
            transmission_efficiency = self.focal_plane_filter * self.entire_filter  # these appear to be the same... not sure

            # equation:
            eff_area = self.geometric_area * reflectance * transmission_efficiency * wavelength * self.instrument_contamination * self.ccd

            self.eff_area = eff_area

        return self.eff_area

    def instrument_response(self, compare=False):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature.
            R(\lambda)=A_{eff}(\lambda,t)G(\lambda)

        effective area A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)
        gain of the CCD-camera system, G(\lambda)=(12398/\lambda/3.65)g
        flat field function F(\mathbf{x})

           Wavelength response functions: calculate the amount of flux per wavelength
         R_i(\lambda) which is equivalent to \eta_i as expressed above.


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

        # TODO: fix doc string here! ^^^
        # TODO: WIP!!   work with values / conversions to get the right number out from this function
        gu = u.count / u.photon
        gain_table = {94: 2.128 * gu, 131: 1.523 * gu, 171: 1.168 * gu, 195: 1.024 * gu, 211: 0.946 * gu,
                      304: 0.658 * gu, 335: 0.596 * gu, 1600: 0.125 * gu, 1700: 0.118 * gu}

        if compare:
            # gives Table 2 values for instruement response
            index = self.wavelength_centers.index(self.wavelength)
            self.inst_response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027, 0.0024, 0.0046][index]
        else:
            # calculations of G from Boerner Table 12
            for key, gain in gain_table.iteritems():
                if str(key) == str(self.wavelength):
                    eff_area = self.effective_area()
                    inst_response = eff_area * gain
                    self.inst_response = inst_response[self.wavelength]

        return self.inst_response

        # NOTES:      ^^^ get the same values as the table when using all table values: self.effective_area(compare_table=True)
        #                 get one of the same values (171) when calculating it with center wavelength - fix


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
