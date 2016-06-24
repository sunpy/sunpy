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

    def __init__(self, **kwargs):
        # no keywords for now. emissivity, temperature, density possible later

        # parameter list
        self.properties = ['wave', 'effarea', 'units', 'geoarea', 'scale', 'platescale',
                           'numfilters', 'wavemin', 'wavestep', 'wavenumsteps', 'wavelog',
                           'filtersincludemesh', 'meshtrans', 'usecontam', 'contamthick',
                           'usephottoelec', 'elecperev', 'usephottodn', 'elecperdn', 'useerror',
                           'fp_filter', 'ent_filter', 'primary', 'secondary', 'ccd', 'contam',
                           'cross_area']

        # wavelength center list
        # TODO: want 6 EUV channels and 2 UV channels --- needs updated
        self.wavelength_centers = [94, 131, 171, 193, 211, 304, 335]

        #           Notes: I needed a way to define self.dataframe locally inside the aia_inst_genx_to_dataframe
        #           and afterwards use it globally. This is what i have found, but am open to other suggestions.

        #           I could import os and recursivly search the .genx directory for the files needed rather than the user defining one file. Then, the variable would change to path_to_genx_dir


        # keyword arguements
        path_to_genx_file = kwargs.get('path_to_genx_file', '')
        save_outfile = kwargs.get('save_outfile', False)
        csv = kwargs.get('csv', '')


        # Load in the filename with the call of the function to initialize channel properties
        if path_to_genx_file and save_outfile:
            self.aia_inst_genx_to_dataframe(path_to_genx_file, save_outfile)
        elif path_to_genx_file:
            self.aia_inst_genx_to_dataframe(path_to_genx_file)
        elif csv:
            self.dataframe = np.array(pd.DataFrame.from_csv(csv))
        else:
            print('Implement keyword path_to_genx_file or csv to load properties.')

        wavelength = kwargs.get('wavelength', '')



        #  These are just for quick plotting for now...
        # assert type(wavelength) == str
        # wavelength range is the same for all channels, but can be loaded per channel, reflectance is only defined if wavelength is defined
        if wavelength:
            if wavelength == 'all':
                wrange = []
                reflect = []
                transm = []
                fp = []     # want to iterate over all properties
                for i in self.wavelength_centers:
                    fp.append(self.dataframe[str(i)]['fp_filter'])
                    wrange.append(self.dataframe[str(i)]['wave'])
                    reflect.append(self.dataframe[str(i)]['primary'] * self.dataframe[str(i)]['secondary'])
                    transm.append(self.dataframe[str(i)]['fp_filter'] * np.array(self.dataframe[str(i)]['ent_filter'])) # * or +?
                self.wavelength_range = wrange
                self.reflectance = reflect
                self.transmittance = transm

            else:
                # add units
                self.wavelength_range = self.dataframe[wavelength]['wave']
                self.reflectance = (self.dataframe[wavelength]['primary']) * (
                    self.dataframe[str(wavelength)]['secondary'])
        else:
            self.wavelength_range = self.dataframe['94']['wave']


            # Notes: I'm trying to define self.stuffs as cleanly as possible - Any input appreciated as this is still new to me
            #          I saw that def __repr__(self): can be used to define self.things but I'm not sure on how..?

    def aia_inst_genx_to_dataframe(self, path_to_file, save=False):
        """
        This definition reads the instrument file aia_V6_all_fullinst, which was obtained from ssw_aia_response_data inside
        ssw_aia_response_genx.tar.gz (an output file saved from SolarSoft Ware (SSW)). It will extract the instrument data and save
        it a dataframe file for easier access.

         Parameters
        ----------
        path_to_file : string, the path location that leads to ssw_aia_response_data/aia_V6_all_fullinst.

        :param save: bool, if True, will save a dataframe containing all of the extracted instrument data to

        Returns
        -------
        :returns self.dataframe: dataframe, stores the dictionary for easy access. Each column lists the channel information and each row is a different peroperty from aia_inst_genx.

        :returns 'channel_properties.csv', dataframe with wavelength centers for each channel is a column and listed in the column are the properties from the is returned where it has keys of the properties of the channel

        Notes:
        np.recarray store information with shape(1,0) and are quite nested.

        """

        # access np.recarray from .genx file
        ssw_array = readsav(path_to_file)
        data = ssw_array['data']

        # store in dataframe
        df = pd.DataFrame()

        # pick out instrument files inside np.recarray
        for name in data.dtype.names:

            # eg.(A94_THICK_FULL or A94_FULL or A94_THICK_FILE)
            # TODO: Fix to read in more than just this file to get all desired channels
            if name.startswith('A') and name.endswith('THICK_FULL'):

                # target number in filename that matches wavelength_center
                start = name.find('A')
                end = name.find('_T')
                wavelength = name[start + 1:end]

                # match the properties in the rec.array to desired properties
                variables = set(data[name][0].dtype.fields.keys()) & set(self.properties)
                assert len(variables) == len(self.properties)

                # then iterate through the properties to fill a dataframe with rec.array information
                for property in variables:
                    try:
                        df.loc[property, wavelength] = data[name][0][property][0]

                    # contam and crossarea are sequences so these don't have to be unpacked as much
                    except ValueError:
                        df.loc[property, wavelength] = data[name][0][property]

        assert len(df) != 0, 'Data Frame is not loading from file.'
        self.dataframe = df

        if save:
            #  save to .csv outfile
            df.to_csv('channel_properties.csv')
            print('saved to channel_properties.csv')

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
        elecperev: float
            a conversion factor for electrons per electronic volts
        ent_filter: array
            transmission efficiency of the entrance filter
        fp_filter: array of floats
            transmission efficiency of the focal plane filter
        contam: array of floats
            correction to adjust for the time varying contamination of the components
        cross_area: array of flats
            cross area of the instrument
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

        # target wavelength specifid and pull out respective properties in dataframe and apply units
        for wavelength, values in self.dataframe.iteritems():
            if int(wavelength) == int(wavelength_center):
                dict['wave'] = values['wave'] * u.angstrom
                dict['effarea'] = values['effarea'] * u.cm ** 2
                dict['units'] = values['units']
                dict['geoarea'] = values['geoarea'] * u.cm ** 2
                dict['scale'] = values['scale']
                dict['platescale'] = values['platescale'] * u.cm ** 2
                dict['numfilters'] = values['numfilters']
                dict['wavemin'] = values['wavemin'] * u.angstrom
                dict['wavestep'] = values['wavestep'] * u.angstrom
                dict['wavenumsteps'] = values['wavenumsteps'] * u.angstrom
                dict['wavelog'] = values['wavelog'] * u.angstrom
                dict['filtersincludemesh'] = values['filtersincludemesh']
                dict['meshtrans'] = values['meshtrans']
                dict['usecontam'] = values['usecontam']
                dict['contamthick'] = values['contamthick']
                dict['usephottoelec'] = values['usephottoelec'] * u.photon/u.electron
                dict['elecperev'] = values['elecperev'] * (u.electron /u.ct) # u.eV?
                dict['usephottodn'] = values['usephottodn'] * u.photon / u.ct
                dict['elecperdn'] = values['elecperdn'] * (u.electron/u.ct)  # the ccd gain!!
                dict['useerror'] = values['useerror']
                dict['fp_filter'] = values['fp_filter']
                dict['ent_filter'] = values['ent_filter']
                dict['primary'] = values['primary']
                dict['secondary'] = values['secondary']
                dict['ccd'] = values['ccd']

                # indexing the next two series did not return expected results because not all have len 6
                index = self.wavelength_centers.index(int(wavelength))
                dict['contam'] = values['contam'][0] * u.cm**2
                dict['cross_area'] = values['cross_area'][index]


                assert len(dict) != 0
                return dict

    def effective_area(self, wavelength,  compare=False):
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

        assert wavelength in self.wavelength_centers, "Choose one integer value from wavelength centers: ' + str(self.wavelength_centers) or 'all'"

        # define variables to call from channel dictionary
        var = self.property_per_channel_dict(wavelength)
        Reflectance = var['primary'] * var['secondary']
        Transmission_efficiency = var['fp_filter'] * var['ent_filter']


        eff_area = var['geoarea'] * Reflectance * Transmission_efficiency * var['contam'] * var['ccd']

        if compare:
            self.eff_area = var['effarea']

        else:
            self.eff_area = eff_area            # print('calc:',wavelength, np.mean(eff_area))
        return self.eff_area

        # TODO: Work on units!! not right as of yet ?!
        # REFERENCE IDL CODE:
        # if usephottoelec and usephottodn:-- not sure why?! - conversions?
        #       units = 12398. / wave * elecperev / elecperdn
        # else:
        #       units = 1
        # eff_area = ((geoarea * ent_filter * fp_filter * primary * secondary * ccd * contam) + cross_area) * units


        #    paper: 0.312 1.172 2.881 1.188 1.206 0.063 0.045


    def instrument_response(self, wavelength = 94, ssw_gain=False, compare = False):
        """

        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature
        \
            R(\lambda)=A_{eff}(\lambda,t)G(\lambda)

        notes:
        For a given position in the image plane \mathbf{x} and wavelength channel i , the pixel values can be expressed as,

        p_i(\mathbf{x})=\int_0^{\infty}\mathrm{d}\lambda\,\eta_i(\lambda)\int_{pixel\,\mathbf{x}}\mathrm{d}\theta\,I(\lambda,\theta)

        Here, \eta_i(\lambda,t,\mathbf{x}) is the efficiency function of the i^{th} channel


        effective area A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)
        gain of the CCD-camera system, G(\lambda)=(12398/\lambda/3.65)g
        flat field function F(\mathbf{x})

           Wavelength response functions: calculate the amount of flux per wavelength
         R_i(\lambda) which is equivalent to \eta_i as expressed above.


            # ccd == quantum efficience ==q
            # calulate gain across wavelength centers - included in master gain of aia image metadata load gain
            # contam = fluctions of the instruement
        Parameters
        ----------
        :keyword
        input: string, area from effective area calculation

        input: may need a photon-to-DN unit conversion   # digital number = counts on detectore

        output: outfile of instrument response per channel

        Returns
        -------
        :return: float, array describing the response per wavelength of effective area (wavelength response)
                want units cm^2 DN phot^-1
        """

        # TODO: fix doc string here! ^^^


        var = self.property_per_channel_dict(wavelength)

        if compare:
            self.inst_response = [0.664, 1.785, 3.365, 1.217, 1.14, 0.041, 0.027]
        else:
            # Calculations of G from Boerner formulas * work in progress*
            if ssw_gain:
                for key, gain in ccd_gain.iteritems():
                    if str(key) == str(wavelength):
                        eff_area = self.effective_area(key)
                        Gain = 12398.0 / (var['wave'] * var['elecperev']) / var['elecperdn'] # < gain
                        inst_response = (eff_area+ var['cross_area']) * Gain ** (u.electron / u.ct)
                        # want: dn/photon  -- gain should be elec/dn??
                        self.inst_response = instr_response

            # calculations of G from Boerner Table 12
            else:
                Gain_dict = {94: 2.128, 131: 1.523, 171: 1.168, 195: 1.024, 211: 0.946, 304: 0.658, 335: 0.596}
                for key, gain in Gain_dict.iteritems():
                    if str(key) == str(wavelength):
                        eff_area = self.effective_area(key)
                        print('eff_area:',eff_area)
                        inst_response = (eff_area  * gain * ((u.electron/u.ct)/var['usephottodn']) * (u.electron / u.ct))# eff_area + var['cross_area'])
                        # print(np.mean(instr_response))
                        self.inst_response = inst_response




        # TODO: get this array to be as expected  , want units: cm**2 DN /phot
        return self.inst_response

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
