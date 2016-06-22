from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__author__ = "Tessa D. Wilkinson"

"""
AIA response functions by integrating ChiantiPy

    - analyze AIA instrument response that tells efficiency of channel using instrument properties
    - then output a wavelength response
    - use spectral model (synthetic spectra) - to get ion emissivities as a function of temp/density
    - obtain temperature response based on chianti spectral contribution functions (emissivity)


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
from astropy import units as u
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
    wave: array
        wavelengths imaged by the instrument and centered around the channel wavelength
    geoarea: float
        the geometric area of the each channel: EUV = 83.0 cm**2, UV = 30.8 cm**2
    elecperev: float
        a conversion factor for electrons per electronic volts
    ent_filter: array
        size of the entire filter on the instrument
    fp_filter: array of floats
        size of fp filter <<<<
    contam: array of floats
        maps the fluctuation of the instrument
    cross_area: array of flats
        cross area of the instrument
    secondary: array of floats
        size of secondary mirror
    primary: array
        size of the primary mirror
    ccd: array, floats
        size of the ccd


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

        # Load in the filename with the call of the function to intialize channel properties
        #    Notes: I needed a way to define self.dataframe locally inside the aia_inst_genx_to_dataframe
        #           and afterwards use it globally. This is what i have found, but am open to other suggestions.

        #   # I can also implement a recursive search of the genx directory if wanted. Then: path_to_genx_dir
        path_to_genx_file = kwargs.get('path_to_genx_file', '')
        save_outfile = kwargs.get('save_outfile', True)
        csv = kwargs.get('csv', '')

        if path_to_genx_file and save_outfile:
            self.aia_inst_genx_to_dataframe(path_to_genx_file, save_outfile)
        elif path_to_genx_file:
            self.aia_inst_genx_to_dataframe(path_to_genx_file)
            print('yeah')
        elif csv:
            self.dataframe = pd.DataFrame.from_csv(csv)
        else:
            print('Implement keyword path_to_genx_file or csv to load properties.')


            # Notes: I'm trying to define self.stuffs as cleanly as possible - Any input appreciated
            #          I saw that def __repr__(self): can be used to define self.things maybe?
            # TODO: research *args **kwargs... more

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

    def properties_per_channel_dictionary(self, wavelength_center=94):
        """

        :param wavelength_center:
        :return:
        """
        assert wavelength_center in self.wavelength_centers, 'Choose one integer value from wavelength centers: ' + str(
            self.wavelength_centers)

        dict = {}

        # target wavelength specifid and pull out respective properties in dataframe
        for wavelength, values in self.dataframe.iteritems():
            if int(wavelength) == int(wavelength_center):
                dict['wave'] = values['wave'] * u.angstrom
                dict['effarea'] = values['effarea'] * u.meter
                dict['units'] = values['units']
                dict['geoarea'] = values['geoarea'] * u.meter
                dict['scale'] = values['scale'] * u.meter
                dict['platescale'] = values['platescale'] * u.meter
                dict['numfilters'] = values['numfilters']
                dict['wavemin'] = values['wavemin'] * u.angstrom
                dict['wavestep'] = values['wavestep'] * u.angstrom
                dict['wavenumsteps'] = values['wavenumsteps'] * u.angstrom
                dict['wavelog'] = values['wavelog'] * u.angstrom
                dict['filtersincludemesh'] = values['filtersincludemesh']
                dict['meshtrans'] = values['meshtrans']
                dict['usecontam'] = values['usecontam']
                dict['contamthick'] = values['contamthick']
                dict['usephottoelec'] = values['usephottoelec'] * u.meter
                dict['elecperev'] = values['elecperev'] * u.meter
                dict['usephottodn'] = values['usephottodn'] * u.meter
                dict['elecperdn'] = values['elecperdn'] * u.meter
                dict['useerror'] = values['useerror']
                dict['fp_filter'] = values['fp_filter'] * u.meter
                dict['ent_filter'] = values['ent_filter'] * u.meter
                dict['primary'] = values['primary'] * u.meter
                dict['secondary'] = values['secondary'] * u.meter
                dict['ccd'] = values['ccd'] * u.meter
                dict['contam'] = values['contam']
                dict['cross_area'] = values['cross_area']

                return dict

    def effective_area(self, wavelength_center):
        """
        Finds the area of the instrument
        needs to work with channels
        AIA instrument response / effective area
            A_{eff}=A_{geo}R_p(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda,t)Q(\lambda)

        Parameters
        ----------
        filename: string, file (or path to file?) with instrument information

        input: a data file giving the area of the instrument

        channel: tell which channel to find the area for, otherwise will just find all channels.

        Returns
        -------
        effective_area: dictionary or array

        """

        if wavelength_center == 'all':
            pass   # run through all
        else:
            assert wavelength_center in self.wavelength_centers, "Choose one integer value from wavelength centers: ' + str(self.wavelength_centers) or 'all'"

            #
            # # replicating the IDL version for now.     idl: if statment usephottoelec and usephottodn-- not sure why?!
            #
            # ones = np.ones(len(wave))
            # print
            # type(wave), type(elecperdn), type(elecperev)
            # units = 12398. / wave * elecperev / elecperdn
            #
            # eff_area = ((geoarea * ent_filter * fp_filter * primary * secondary * ccd * contam) + cross_area) * units
            #
            # return eff_area

    def instrument_response(self, eff_area, gain, flat_field):
        """
        Describes the instrument (AIA) response per wavelength by calculating effective area vs wavelength of the strongest emission lines present in the solar feature
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

        """
        # TODO: implement response function and plot!
        pass


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
