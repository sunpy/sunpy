"""
Calculate SDO/AIA wavelength and temperature response functions.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import warnings
from collections import namedtuple

import numpy as np
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import splev, splrep
try:
    from ChiantiPy.tools.data import MasterList
except ImportError:
    MasterList = None

from sunpy.util.config import get_and_create_download_dir
from sunpy.util.net import check_download_file
from .response_utils import aia_instr_properties_to_table, make_emiss_table, EmissTableInterface

SSW_AIA_REMOTE_PATH = 'https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/'
AIA_INSTR_FILES = ['aia_V6_all_fullinst.genx', 'aia_V6_fuv_fullinst.genx']

__author__ = ["Tessa D. Wilkinson", "Will Barnes"]


class Response(object):
    """
    Calculate the SDO/AIA wavelength and temperature response functions.

    This class calculates the wavelength and temperature response functions for
    the channels of the Atmospheric Imaging Assembly (AIA) instrument on
    the Solar Dynamics Observatory spacecraft. The response functions are
    calculated using AIA instrument information obtained from
    Solar Software (SSW) .genx files. Temperature response functions are
    calculated using atomic data from the CHIANTI atomic database. See [1]_ for
    more details.

    Parameters
    ----------
    channel_list : `list`, optional
        AIA channels, defaults to the 7 EUV channels, [94,]
    instrument_files : `list`, optional
        If no instrument files are supplied and they do not exist locally, the
        necessary files are automatically downloaded.

    Examples
    --------
    >>> from sunpy.instr.aia import Response
    >>> response = Response()
    >>> wavelength_response = response.calculate_wavelength_response()
    >>> temperature_response = response.calculate_temperature_response()

    Notes
    -----
    Calculating the temperature response functions requires tabulated emissivities,
    for both the continuum and line contributions, for each ion. The first time the temperature
    response functions are calculated, this table is automatically generated. You can also
    generate your own emissivity table using `make_emiss_table` and then passing the appropriate
    filename to the temperature response function method.

    References
    ----------
    .. [1] Boerner et al., 2012, Sol. Phys., `275, 41
        <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`_
    """

    def __init__(self, channel_list=None, instrument_files=None):
        if channel_list is None:
            channel_list = [94, 131, 171, 193, 335, 211, 304]
        self._get_channel_info(channel_list, instrument_files)

    def _get_channel_info(self, channel_list, instrument_files):
        """
        Get instrument info for channels in channel list.

        Download the instrument files if needed or if they were not specified at the command
        line. Creates the dictionary `self._channel_info`.
        """
        if instrument_files is None:
            instrument_files = []
            download_dir = get_and_create_download_dir()
            for instr_file in AIA_INSTR_FILES:
                check_download_file(instr_file, SSW_AIA_REMOTE_PATH, download_dir=download_dir)
                instrument_files.append(os.path.join(download_dir, instr_file))

        data_table = aia_instr_properties_to_table(channel_list, instrument_files)
        self._channel_info = {}
        for c in channel_list:
            index = channel_list.index(c)
            self._channel_info[c] = {property: data_table[property][index]
                                     for property in data_table.columns.keys()}

    def calculate_effective_area(self, channel, include_crosstalk=True, **kwargs):
        """
        Calculate the effective area for a given channel of the instrument.

        According to [1]_, the effective area is given by,

        .. math::
            A_{eff}(\lambda,t) = A_{geo}R_P(\lambda)R_S(\lambda)T_E(\lambda)
            T_F(\lambda)D(\lambda,t)Q(\lambda),

        in units of :math:`\mathrm{cm}^2`, where,

        - :math:`A_{geo}`: geometrical collecting area
        - :math:`R_P`, :math:`R_S`: reflectances of primary and secondary mirrors, respectively
        - :math:`T_E`, :math:`T_F`: transmission efficiency of the entrance and focal-plane filters, respectively
        - :math:`D`: correction to account for time-dependent contamination and deterioration
        - :math:`Q`: quantum efficiency of the CCD

        The effective area contains information about the efficiency of the
        telescope optics and its sensitivity as a function of wavelength. All of the telescope
        properties are read from the AIA instrument files available in SolarSoft.

        Parameters
        ----------
        channel : `int`
            wavelength center of the channel
        include_crosstalk : `bool`

        Returns
        -------
        effective_area : array-like
            effective area of the instrument over the entire wavelength range
        """
        reflectance = (self._channel_info[channel]['primary_mirror_reflectance']
                       * self._channel_info[channel]['secondary_mirror_reflectance'])
        transmission_efficiency = (self._channel_info[channel]['focal_plane_filter_efficiency']
                                   * self._channel_info[channel]['entrance_filter_efficiency'])
        effective_area = (self._channel_info[channel]['geometric_area_ccd']
                          * reflectance*transmission_efficiency
                          * self._channel_info[channel]['ccd_contamination']
                          * self._channel_info[channel]['quantum_efficiency_ccd'])
        if include_crosstalk:
            effective_area += self.calculate_crosstalk(channel)

        return effective_area

    def calculate_crosstalk(self, channel):
        """
        Calculate effects due to channel crosstalk.

        The focal-plane filters on Telescopes 1 (131 and 335) and 4 (94 and 304)
        do not perfectly reject light from the opposite channels. This function
        implements a correction for the given channel, based on the method in
        SSW, specifically that in `sdo/aia/idl/response/aia_bp_blend_channels.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_blend_channels.pro>`_.
        See section 2.2.1 and Figure 10 of [1]_ for more details.
        """
        # check which channel is the contaminator
        if channel == 94:
            cross_channel = 304
        elif channel == 131:
            cross_channel = 335
        elif channel == 304:
            cross_channel = 94
        elif channel == 335:
            cross_channel = 131
        else:
            cross_channel = None

        # check to see if the contamination can/should be calculated
        if cross_channel not in self._channel_info or cross_channel is None:
            if cross_channel is not None:
                warnings.warn(('{cross} is not included in the channel list. '
                               'Cross-contamination between {cross} and '
                               '{primary} will not be included.').format(
                                cross=cross_channel, primary=channel))
            return (np.zeros(len(self._channel_info[channel]['wavelength']))
                    * self._channel_info[channel]['effective_area'].unit)

        return (self._channel_info[cross_channel]['primary_mirror_reflectance']
                * self._channel_info[cross_channel]['secondary_mirror_reflectance']
                * self._channel_info[cross_channel]['quantum_efficiency_ccd']
                * self._channel_info[cross_channel]['ccd_contamination']
                * self._channel_info[cross_channel]['geometric_area_ccd']
                * self._channel_info[cross_channel]['entrance_filter_efficiency']
                * self._channel_info[channel]['focal_plane_filter_efficiency'])

    def calculate_system_gain(self, channel):
        """
        Calculate the gain of the CCD camera system for a given channel in units of DN
        :math:`\mathrm{photon}^{-1}`.

        Parameters
        ----------
        channel : `int`
            the wavelength center of the channel
        """
        # convert hc to convenient units
        hc = (const.h * const.c).to(u.eV * u.angstrom)
        # energy per photon in units of eV
        ev_per_photon = hc/self._channel_info[channel]['wavelength']/u.photon
        electron_per_ev = self._channel_info[channel]['electron_per_ev']
        electron_per_photon = ev_per_photon*electron_per_ev
        gain = electron_per_photon/self._channel_info[channel]['electron_per_dn']
        return gain

    def calculate_wavelength_response(self, include_crosstalk=True):
        """
        Calculate the wavelength response function for all channels of the
        instrument.

        The wavelength response function for channel :math:`i` is the product of
        the effective area and the instrument gain,

        .. math::
            R_i(\lambda) = A_{eff}(\lambda)G(\lambda)

        in units of :math:`\mathrm{cm}^2` DN :math:`\mathrm{photon}^{-1}`.

        Parameters
        ----------
        include_crosstalk : `bool`
            See `calculate_crosstalk`
        """
        Response = namedtuple('Response', 'wavelength response')
        wavelength_response = {}
        for channel in self._channel_info:
            system_gain = self.calculate_system_gain(channel)
            effective_area = self.calculate_effective_area(channel,
                                                           include_crosstalk=include_crosstalk)
            wavelength_response[channel] = Response(
                wavelength=self._channel_info[channel]['wavelength'],
                response=system_gain*effective_area)

        return wavelength_response

    def calculate_temperature_response(self, ion_list=None, emiss_table_file=None,
                                       include_crosstalk=True, **kwargs):
        """
        Calculate the temperature response functions.

        According to [1]_, the temperature response function for channel :math:`i` is
        given by,

        .. math::
           K_i(T) = \int_0^{\infty}\mathrm{d}\lambda\,G(\lambda,n,T)R_i(\lambda),

        in units of DN :math:`\mathrm{cm}^{5}\,\mathrm{s}^{-1}\,\mathrm{pixel}^{-1}` where
        :math:`R_i(\lambda)` is the wavelength response function of channel :math:`i` and
        :math:`G(\lambda,T)` is the contribution function. Note that :math:`G` includes both
        the line and continuum contributions such that,

        .. math::
           G(\lambda,n,T) = G_{ff}(\lambda,T) + G_{fb}(\lambda,T) + G_{tp}(\lambda,T) + \sum_{Z,z}\\frac{1}{4\pi}\mathrm{Ab}(Z)\\frac{N(Z^{+z})}{N(Z)}\\varepsilon_{Z,z}(\lambda,n,T)\\frac{1}{n}

        where :math:`G_{ff}(\lambda,T),\,G_{fb}(\lambda,T),\,G_{tp}(\lambda,T)` are the free-free,
        free-bound, and two-photon continuum emissivities.

        Parameters
        ----------
        ion_list : array-like, optional
            List of ion names, e.g. 'fe_15' for Fe XV. Defaults to all ions in CHIANTI database
        emiss_table_file : `str`, optional
            Defaults to 'aia_emiss_table.h5' in the default download directory for SunPy


        .. warning::
            The first time you calculate the temperature response functions, you will
            need to build the emission table by calculating the contribution function for each
            ion. This may take up to 30 minutes, but will only need to be done once.
        """
        if emiss_table_file is None:
            emiss_table_file = os.path.join(get_and_create_download_dir(), 'aia_emiss_table.h5')
        if not os.path.exists(emiss_table_file):
            warnings.warn('Building emissivity table {}. This may take a while, but only needs to be done once.'.format(emiss_table_file))
            make_emiss_table(emiss_table_file, **kwargs)

        wavelength_response = self.calculate_wavelength_response(include_crosstalk=include_crosstalk)
        table_interface = EmissTableInterface(emiss_table_file)
        tmp_tresponse = {channel: np.zeros(table_interface.temperature.shape) for channel in wavelength_response}

        if ion_list is None:
            if MasterList is None:
                raise ImportError('CHIANTI ion list not available. Install ChiantiPy or provide a list of ions.')
            ion_list = MasterList

        for ion in ion_list:
            tmp_ion = table_interface[ion]
            for channel in wavelength_response:
                # interpolate wavelength response
                wave_interp = splrep(wavelength_response[channel].wavelength.value,
                                     wavelength_response[channel].response.value)
                # line emission
                rsp_func = splev(tmp_ion.transitions.value, wave_interp, ext=1)
                rsp_func = np.where(rsp_func < 0.0, 0.0, rsp_func)
                line_contrib = (np.dot(tmp_ion.contribution_function.value, rsp_func)
                                * tmp_ion.contribution_function.unit
                                * wavelength_response[channel].response.unit)
                # continuum
                cont_total = tmp_ion.two_photon_continuum + tmp_ion.free_free_continuum + tmp_ion.free_bound_continuum
                rsp_func = splev(table_interface.continuum_wavelength.value, wave_interp)
                rps_func = np.where(rsp_func < 0.0, 0.0, rsp_func)*np.gradient(table_interface.continuum_wavelength.value)
                cont_contrib = (np.dot(cont_total.value, rsp_func)
                                * cont_total.unit*table_interface.continuum_wavelength.unit
                                * wavelength_response[channel].response.unit)
                # total
                if not hasattr(tmp_tresponse[channel],'unit'):
                    tmp_tresponse[channel] = (tmp_tresponse[channel]*line_contrib.unit
                                              * self._channel_info[channel]['plate_scale'].unit)
                tmp_tresponse[channel] += (line_contrib + cont_contrib)*self._channel_info[channel]['plate_scale']

        Response = namedtuple('Response', 'temperature density response')
        return {channel: Response(temperature=table_interface.temperature,
                                  density=table_interface.density,
                                  response=tmp_tresponse[channel])
                for channel in tmp_tresponse}


