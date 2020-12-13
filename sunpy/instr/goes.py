"""
Contains functions useful for analysing GOES/XRS data.

Each of the Geostationary Operational Environmental Satellite (GOES) series
since the mid-1970s has carried an X-Ray Sensor (XRS) which observes
full-disk-integrated solar flux in two broadband channels:
1--8 angstrom (long); and 0.5--4 angstrom (short). For more information on
the GOES/XRS instrument, see [Ref1]_.  GOES/XRS has become
the "standard candle" for solar coronal observations due its longevity and
consistency. The GOES event list, based on GOES/XRS observations, has
become the standard solar flare catalogue.
See https://www.ngdc.noaa.gov/stp/solar/solarflares.html for information
on the GOES event list definitions and data.

The functions in this module provide useful software to analyse GOES/XRS
observations. First they allow the GOES event list to be imported into a
Python session (`~sunpy.instr.goes.get_goes_event_list`).

They also allow the thermodynamic properties of the emitting solar plasma to be
determined. Temperature and emission measure are obtained using
`~sunpy.instr.goes.calculate_temperature_em`, this function currently relies on
lookup tables relating the GOES fluxes to the isothermal temperature and volume
emission measure. These tables were calculated by functions in SolarSoftWare
(SSW) using the CHIANTI atomic physics database ([Ref2]_). For more detail, see
the docstring of calculate_temperature_em` and references therein.

The radiative loss rate of the soft X-ray-emitting plasma across all
wavelengths can be found with
`~sunpy.instr.goes.calculate_radiative_loss_rate`, which makes use of a look up
table calculated by functions in SSW using CHIANTI. This table relates the
temperature and emission measure of the emitting solar plasma to the thermal
energy radiated over all wavelengths. For more information on how this is
done, see the docstring of `~sunpy.instr.goes._calc_rad_loss` and reference
therein.

Meanwhile, the X-ray luminosity in the two GOES passbands can be obtained by
`~sunpy.instr.goes.calculate_xray_luminosity`. To do so, this function calls
`~sunpy.instr.goes._goes_lx` and `~sunpy.instr.goes.calc_xraylum`.

References
----------

.. [Ref1] Hanser, F.A., & Sellers, F.B. 1996, Proc. SPIE, 2812, 344
.. [Ref2] Dere, K.P., et al. 2009 A&A, 498, 915 DOI:
          `10.1051/0004-6361/200911712 <https://doi.org/10.1051/0004-6361/200911712>`__
"""

import csv
import copy
import socket
import datetime
from itertools import dropwhile
from urllib.parse import urljoin

import numpy as np
from scipy import interpolate
from scipy.integrate import cumtrapz, trapz

import astropy.units as u
from astropy.time import TimeDelta

from sunpy import timeseries
from sunpy.coordinates import sun
from sunpy.data import manager
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.util.config import get_and_create_download_dir

GOES_CONVERSION_DICT = {'X': u.Quantity(1e-4, "W/m^2"),
                        'M': u.Quantity(1e-5, "W/m^2"),
                        'C': u.Quantity(1e-6, "W/m^2"),
                        'B': u.Quantity(1e-7, "W/m^2"),
                        'A': u.Quantity(1e-8, "W/m^2")}
__all__ = ['get_goes_event_list', 'calculate_temperature_em',
           'calculate_radiative_loss_rate', 'calculate_xray_luminosity', 'flux_to_flareclass',
           'flareclass_to_flux']

try:
    # Check required data files are present in user's default download dir
    # Define location where GOES data files are stored.
    # Manually resolve the hostname
    HOST = socket.gethostbyname_ex('hesperia.gsfc.nasa.gov')[0]
except socket.gaierror:
    HOST = ''
GOES_REMOTE_PATH = f"http://{HOST}/ssw/gen/idl/synoptic/goes/"
# Define variables for file names
FILE_TEMP_COR = "goes_chianti_temp_cor.csv"
FILE_TEMP_PHO = "goes_chianti_temp_pho.csv"
FILE_EM_COR = "goes_chianti_em_cor.csv"
FILE_EM_PHO = "goes_chianti_em_pho.csv"
FILE_RAD_COR = "chianti7p1_rad_loss.txt"


def get_goes_event_list(timerange, goes_class_filter=None):
    """
    Retrieve list of flares detected by GOES within a given time range.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        The time range to download the event list for.
    goes_class_filter: `str`, optional
        A string specifying a minimum GOES class for inclusion in the list,
        e.g., "M1", "X2".

    Returns
    -------
    `list`:
        A list of all the flares found for the given time range.
    """
    # Importing hek and attrs here to avoid calling code that relies on optional dependencies.
    from sunpy.net import attrs as a
    from sunpy.net import hek

    # use HEK module to search for GOES events
    client = hek.HEKClient()
    event_type = 'FL'
    tstart = timerange.start
    tend = timerange.end

    # query the HEK for a list of events detected by the GOES instrument
    # between tstart and tend (using a GOES-class filter)
    if goes_class_filter:
        result = client.search(a.Time(tstart, tend),
                               a.hek.EventType(event_type),
                               a.hek.FL.GOESCls > goes_class_filter,
                               a.hek.OBS.Observatory == 'GOES')
    else:
        result = client.search(a.Time(tstart, tend),
                               a.hek.EventType(event_type),
                               a.hek.OBS.Observatory == 'GOES')

    # want to condense the results of the query into a more manageable
    # dictionary
    # keep event data, start time, peak time, end time, GOES-class,
    # location, active region source (as per GOES list standard)
    # make this into a list of dictionaries
    goes_event_list = []

    for r in result:
        goes_event = {
            'event_date': parse_time(r['event_starttime']).strftime(
                '%Y-%m-%d'),
            'start_time': parse_time(r['event_starttime']),
            'peak_time': parse_time(r['event_peaktime']),
            'end_time': parse_time(r['event_endtime']),
            'goes_class': str(r['fl_goescls']),
            'goes_location': (r['event_coord1'], r['event_coord2']),
            'noaa_active_region': r['ar_noaanum']
        }
        goes_event_list.append(goes_event)

    return goes_event_list


def calculate_temperature_em(goests, abundances="coronal",
                             download=False, download_dir=None):
    """
    Calculates temperature and emission measure from a
    `~sunpy.timeseries.sources.XRSTimeSeries`.

    This function calculates the isothermal temperature and
    corresponding volume emission measure of the solar soft X-ray
    emitting plasma observed by the GOES/XRS. This is done using the
    observed flux ratio of the short (0.5-4 angstrom) to long (1-8 angstrom)
    channels. The results are returned in a new `~sunpy.timeseries.TimeSeries` object which
    contains metadata and flux data of the input TimeSeries object in
    addition to the newly found temperature and emission measure values.

    Parameters
    ----------
    goeslc : `~sunpy.timeseries.sources.XRSTimeSeries`
        The TimeSeries containing GOES flux data which **MUST**
        be in units of "W/m^2".
    abundances : {'coronal' | 'photospheric'}, optional
        States whether "photospheric" or "coronal" abundances should be
        assumed, default to 'coronal'.
    download : `bool`, optional
        If `True`, the GOES temperature and emission measure data files are
        downloaded. It is important to do this if a new version of the files
        has been generated due to a new CHIANTI version being released or the
        launch of new GOES satellites since these files were last downloaded.
        Defaults to `False`.
    download_dir : `str`, optional
        The directory to download the GOES temperature and emission measure
        data files to, defaults to the default download directory.

    Returns
    -------
    `~sunpy.timeseries.sources.XRSTimeSeries`
        Contains same metadata and data as input timeseries with the
        following two additional data columns:

        | ts_new.to_dataframe().temperature - Array of temperatures [MK]
        | ts_new.to_dataframe().em - Array of volume emission measures [cm**-3]

    Notes
    -----
    The temperature and volume emission measure are calculated here
    using the methods of White et al. (2005) who used the
    CHIANTI atomic physics database to model the response of the ratio
    of the short (0.5-4 angstrom) to long (1-8 angstrom) channels of the
    XRSs on board various GOES satellites. This method assumes an
    isothermal plasma, the ionization equilibria of
    [2]_, and a constant density of 10**10 cm**-3.
    (See [1]_ for justification of this last assumption.)
    This function is based on "goes_chianti_tem.pro" in SolarSoftWare
    written in IDL by Stephen White.

    Recent fluxes released to the public are scaled to be consistent
    with GOES-7. In fact these recent fluxes are correct and so this
    correction must be removed before proceeding to use transfer
    functions.

    Measurements of short channel flux of less than 1e-10 W/m**2 or
    long channel flux less than 3e-8 W/m**2 are not considered good.
    Ratio values corresponding to such fluxes are set to 0.003.

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231, DOI: 10.1007/s11207-005-2445-z
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> import sunpy.timeseries as ts
    >>> from sunpy.instr.goes import calculate_temperature_em
    >>> from sunpy.data.sample import GOES_XRS_TIMESERIES  # doctest: +REMOTE_DATA
    >>> goests = ts.TimeSeries(GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> goests.to_dataframe()[0:10]  # doctest: +REMOTE_DATA
                                           xrsa          xrsb
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07
    >>> goests_new = calculate_temperature_em(goests)  # doctest: +REMOTE_DATA
    >>> goests_new.to_dataframe()[0:10]  # doctest: +REMOTE_DATA
                                           xrsa          xrsb  temperature            em
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07     3.503510  2.190626e+48
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07     3.534262  2.055847e+48
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07     3.518700  2.122771e+48
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07     3.550100  1.990333e+48
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07     3.518700  2.122771e+48
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07     3.550100  1.990333e+48
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07     3.518700  2.122771e+48
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07     3.534262  2.055847e+48
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07     3.550100  1.990333e+48
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07     3.534262  2.055847e+48
    """
    # Check that input argument is of correct type
    if not isinstance(goests, timeseries.XRSTimeSeries):
        raise TypeError("goests must be a XRSTimeSeries object")
    if not download_dir:
        download_dir = get_and_create_download_dir()

    # Find temperature and emission measure with _goes_chianti_tem
    temp, em = _goes_chianti_tem(
        goests.quantity("xrsb"),
        goests.quantity("xrsa"),
        satellite=goests.meta.metas[0]["TELESCOP"].split()[1],
        date=goests.to_dataframe().index[0],
        abundances=abundances, download=download, download_dir=download_dir)

    ts_new = timeseries.XRSTimeSeries(meta=copy.deepcopy(goests.meta),
                                      data=copy.deepcopy(goests.to_dataframe()),
                                      units=copy.deepcopy(goests.units))
    ts_new = ts_new.add_column("temperature", temp)
    ts_new = ts_new.add_column("em", em)

    return ts_new


@u.quantity_input
def _goes_chianti_tem(longflux: u.W/u.m/u.m, shortflux: u.W/u.m/u.m, satellite=8,
                      date=datetime.datetime.today(), abundances="coronal",
                      download=False, download_dir=None):
    """
    Calculates temperature and emission measure from GOES/XRS data.

    This function calculates the isothermal temperature and volume
    emission measure of the solar soft X-ray emitting plasma observed by
    the GOES/XRS.  This is done using the observed flux ratio of the
    short (0.5-4 angstrom) to long (1-8 angstrom) channels.

    Parameters
    ----------
    longflux, shortflux : `~astropy.units.Quantity`
        Arrays containing the long and short GOES/XRS flux measurements
        respectively as a function of time.  Must be of same length. [W/m**2].

    satellite : int (optional)
        Number of GOES satellite used to make observations, important for
        correct calibration of data.
        Default=8

    date : `astropy.time.Time` or `str`
        Date when observations made.  Important for correctcalibration.
        Default=today

    abundances : (optional) string equalling 'coronal' or 'photospheric'
        States whether photospheric or coronal abundances should be
        assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES temperature and emission measure data files are
        downloaded.  It is important to do this if a new version of the files
        has been generated due to a new CHIANTI version being released or the
        launch of new GOES satellites since these files were last downloaded.
        Default=False

    download_dir : (optional) str
        The directory to download the GOES temperature and emission measure
        data files to.
        Default=SunPy default download directory

    Returns
    -------
    temp : `~astropy.units.Quantity`
        Array of temperature values of same length as longflux and
        shortflux. Units=[MK]

    em : `~astropy.units.Quantity`
        Array of volume emission measure values of same length as longflux
        and shortflux.  Units=[10**49 cm**-3]

    Notes
    -----
    The temperature and volume emission measure are calculated here
    using the methods of [1]_ who used the
    CHIANTI atomic physics database to model the response of the ratio
    of the short (0.5-4 angstrom) to long (1-8 angstrom) channels of the
    XRSs on board various GOES satellites.  This method assumes an
    isothermal plasma, the ionisation equilibria of
    [2]_, and a constant density of 10**10 cm**-3.
    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_chianti_tem.pro in SolarSoftWare
    written in IDL by Stephen White.

    Recent fluxes released to the public are scaled to be consistent
    with GOES-7.  In fact these recent fluxes are correct and so this
    correction must be removed before proceeding to use transfer
    functions.
    Email Rodney Viereck (NOAA) for more information.

    Measurements of short channel flux of less than 1e-10 W/m**2 or
    long channel flux less than 3e-8 W/m**2 are not considered good.
    Ratio values corresponding to such fluxes are set to 0.003.

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231, DOI: 10.1007/s11207-005-2445-z
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> from sunpy.instr.goes import _goes_chianti_tem
    >>> from astropy.units import Quantity
    >>> longflux = Quantity([7e-6, 7e-6], unit="W/m**2")
    >>> shortflux = Quantity([7e-7, 7e-7], unit="W/m**2")
    >>> temp, em = _goes_chianti_tem(longflux, shortflux, satellite=15,
    ...                              date='2014-04-16',
    ...                              abundances="coronal")  # doctest: +REMOTE_DATA
    >>> temp  # doctest: +REMOTE_DATA
    <Quantity [11.28295376, 11.28295376] MK>
    >>> em  # doctest: +REMOTE_DATA
    <Quantity [4.78577516e+48, 4.78577516e+48] 1 / cm3>
    """
    if not download_dir:
        download_dir = get_and_create_download_dir()
    # ENSURE INPUTS ARE OF CORRECT TYPE AND VALID VALUES
    longflux = longflux.to(u.W/u.m/u.m)
    shortflux = shortflux.to(u.W/u.m/u.m)
    satellite = int(satellite)
    if satellite < 1:
        raise ValueError("satellite must be the number of a "
                         "valid GOES satellite (>1).")
    date = parse_time(date)
    # Check flux arrays are of same length.
    if len(longflux) != len(shortflux):
        raise ValueError(
            "longflux and shortflux must have same number of elements.")

    # PREPARE DATA
    # GOES 6 long channel flux before 1983-Jun-28 must be corrected by a
    # factor of 4.43/5.32
    if date < parse_time((1983, 6, 28)) and satellite == 6:
        longflux_corrected = longflux*(4.43/5.32)
    else:
        longflux_corrected = longflux
    # Un-scale fluxes if GOES satellite is after 7.  See 2nd paragraph
    # in Notes section of docstring above.
    if satellite > 7:
        longflux_corrected = longflux_corrected / 0.7
        shortflux_corrected = shortflux / 0.85
    else:
        shortflux_corrected = shortflux
    # Calculate short to long channel ratio.
    # Data which is not good have their ratio value set to 0.003.
    # See Notes section in docstring above.
    index = np.logical_or(
        shortflux_corrected < u.Quantity(1e-10, unit="W/m**2"),
        longflux_corrected < u.Quantity(3e-8, unit="W/m**2"))
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio.value[index] = u.Quantity(0.003, unit="W/m**2")

    # FIND TEMPERATURE AND EMISSION MEASURE FROM FUNCTIONS BELOW
    temp = _goes_get_chianti_temp(fluxratio, satellite=satellite,
                                  abundances=abundances, download=download,
                                  download_dir=download_dir)
    em = _goes_get_chianti_em(longflux_corrected, temp, satellite=satellite,
                              abundances=abundances, download=download,
                              download_dir=download_dir)
    return temp, em


@manager.require('file_temp_cor',
                 [urljoin(GOES_REMOTE_PATH, FILE_TEMP_COR)],
                 '3d8ddaaabf0faf75ba8d15e0c468896ce3d7622cc23076bf91437951e0ab3ad4')
@manager.require('file_temp_pho',
                 [urljoin(GOES_REMOTE_PATH, FILE_TEMP_PHO)],
                 'dd8c6b949a492174146a0b7307dd5fb197236431dbbedfdbab2e3f8dcd360267')
@u.quantity_input
def _goes_get_chianti_temp(fluxratio: u.one, satellite=8, abundances="coronal",
                           download=False, download_dir=None):
    """
    Calculates temperature from GOES flux ratio.

    This function calculates the isothermal temperature of the solar
    soft X-ray emitting plasma observed by the GOES/XRS from the
    observed flux ratio of the short (0.5-4 angstrom) to
    long (1-8 angstrom) channels.  This function is not intended to be
    called directly but by _goes_chianti_tem(), although it can be used
    independently.  However, if used independently data preparation,
    such as correctly rescaling fluxes for some satellites etc. will
    not be carried out.  This is done in _goes_chianti_tem().

    Parameters
    ----------
    fluxratio : `~astropy.units.Quantity`
        Array containing the ratio of short channel to long channel
        GOES/XRS flux measurements.

    satellite : int (optional)
        Number of GOES satellite used to make observations. Important for
        correct calibration of data.
        Default=8

    abundances : (optional) string equalling 'coronal' or 'photospheric'
        States whether photospheric or coronal abundances should be
        assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES temperature data files are downloaded.
        It is important to do this if a new version of the files has been
        generated due to a new CHIANTI version being released or the launch
        of new GOES satellites since these files were last downloaded.
        Default=False

    download_dir : (optional) str
        The directory to download the GOES temperature data file to.
        Default=SunPy default download directory

    Returns
    -------
    temp : `~astropy.units.Quantity`
        Array of temperature values of same length as longflux and
        shortflux. Units=[MK]

    Notes
    -----
    This function uses csv files representing the modelled relationship
    between temperature of the soft X-ray emitting plasma and the
    short to long channel GOES flux ratio.  goes_chianti_temp_cor.csv
    is used when coronal abundances are assumed while
    goes_chianti_temp_pho.csv is used when photospheric abundances are
    assumed.  (See make_goes_chianti_temp.py for more detail.)

    These files were calculated using the methods of [1]_
    who used the CHIANTI atomic physics database to model the response
    of the ratio of the short (0.5-4 angstrom) to long (1-8 angstrom)
    channels of the XRSs on board various GOES satellites.  This method
    assumes an isothermal plasma, the ionisation equilibria of

    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_get_chianti_temp.pro in
    SolarSoftWare written in IDL by Stephen White.

    For correct preparation of GOES data before calculating temperature
    see _goes_chianti_tem() (Notes section of docstring).

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231, DOI: 10.1007/s11207-005-2445-z
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> from astropy.units import Quantity
    >>> from sunpy.instr.goes import _goes_get_chianti_temp
    >>> fluxratio = Quantity([0.1,0.1])
    >>> temp = _goes_get_chianti_temp(fluxratio, satellite=15,
    ...                               abundances="coronal")  # doctest: +REMOTE_DATA
    >>> temp  # doctest: +REMOTE_DATA
    <Quantity [12.27557778, 12.27557778] MK>
    """
    # check inputs are correct
    fluxratio = fluxratio.decompose()
    int(satellite)
    if satellite < 1:
        raise ValueError("satellite must be the number of a "
                         "valid GOES satellite (>1).")
    # if abundance input is valid create file suffix, abund, equalling
    # of 'cor' or 'pho'.
    if abundances == "coronal":
        data_file = manager.get('file_temp_cor')
    elif abundances == "photospheric":
        data_file = manager.get('file_temp_pho')
    else:
        raise ValueError("abundances must be a string equalling "
                         "'coronal' or 'photospheric'.")

    # Initialize lists to hold model data of flux ratio - temperature
    # relationship read in from csv file
    modeltemp = []  # modelled temperature is in log_10 space in units of MK
    modelratio = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = f"ratioGOES{satellite}"
    # Read data representing appropriate temperature--flux ratio
    # relationship depending on satellite number and assumed abundances.
    with open(data_file, "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelratio.append(float(row[label]))
    modeltemp = np.asarray(modeltemp)
    modelratio = np.asarray(modelratio)

    # Ensure input values of flux ratio are within limits of model table
    if np.min(fluxratio) < np.min(modelratio) or \
       np.max(fluxratio) > np.max(modelratio):
        raise ValueError(
            "For GOES {0}, all values in fluxratio input must be within " +
            "the range {1} - {2}.".format(satellite, np.min(modelratio),
                                          np.max(modelratio)))

    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modelratio, modeltemp, s=0)
    temp = 10.**interpolate.splev(fluxratio.value, spline, der=0)
    temp = u.Quantity(temp, unit='MK')

    return temp


@manager.require('file_em_cor',
                 [urljoin(GOES_REMOTE_PATH, FILE_EM_COR)],
                 'a7440e20cbcb74e87db528e8e9d47cd69fbbd8f56ddc92cf4e854a66fb2a6172')
@manager.require('file_em_pho',
                 [urljoin(GOES_REMOTE_PATH, FILE_EM_PHO)],
                 '0d59042b265bf76351d129b3e2a5844b3a9c96943cb246538013fd8c1b9b71b9')
@u.quantity_input
def _goes_get_chianti_em(longflux: u.W/u.m/u.m, temp: u.MK, satellite=8,
                         abundances="coronal", download=False,
                         download_dir=None):
    """
    Calculates emission measure from GOES 1-8A flux and temperature.

    This function calculates the emission measure of the solar
    soft X-ray emitting plasma observed by the GOES/XRS from the
    the ratio of the isothermal temperature and observed long channel
    (1-8 angstrom) flux which scales with the emission measure.
    This function is not intended to be called directly but by
    _goes_chianti_tem(), although it can be used independently.
    However, if used independently data preparation, such as correctly
    rescaling fluxes for some satellites etc. will not be carried out.
    This is done in _goes_chianti_tem().

    Parameters
    ----------
    longflux : `~astropy.units.Quantity`
        Array containing the observed GOES/XRS long channel flux.
        Units=[W/m**2]

    temp : `~astropy.units.Quantity`
        Array containing the GOES temperature.  Units=[MK]

    satellite : int (optional)
        Number of GOES satellite used to make observations.
        Important for correct calibration of data.
        Default=8

    abundances : (optional) {'coronal' | 'photospheric'}
        States whether photospheric or coronal abundances should be
        assumed.
        Default='coronal'

    download : (optional) `bool`
        If True, the GOES emission measure data file is downloaded.
        It is important to do this if a new version of the file has been
        generated due to a new CHIANTI version being released or the launch of
        new GOES satellites since these file was last downloaded.
        Default=False

    download_dir : (optional) `str`
        The directory to download the GOES emission measure data file to.
        Default=SunPy default download directory

    Returns
    -------
    em : `~astropy.units.Quantity`
         Array of emission measure values of same length as longflux
         and temp.  [cm**-3]

    Notes
    -----
    This function uses csv files representing the modelled relationship
    between the temperature of the solar soft X-ray emitting plasma
    and the resulting observed flux in the GOES/XRS long channel
    (1-8 angstroms).  goes_chianti_em_cor.csv is used when coronal
    abundances are assumed while goes_chianti_em_pho.csv is used when
    photospheric abundances are assumed.
    (See make_goes_chianti_temp.py for more detail.)

    These files were calculated using the methods of White et al. (2005)
    who used the CHIANTI atomic physics database and GOES transfer
    functions to model the response of the long channel to the
    temperature of the emitting plasma for XRSs on board various GOES
    satellites.  The emission measure can then be found by scaling the
    ratio of these two properties.  This method assumes an isothermal
    plasma, the ionisation equilibria of Mazzotta et al. (1998), and
    a constant density of 10**10 cm**-3.
    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_get_chianti_temp.pro in
    SolarSoftWare written in IDL by Stephen White.

    For correct preparation of GOES data before calculating temperature
    see _goes_chianti_tem() (Notes section of docstring).

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231, DOI: 10.1007/s11207-005-2445-z
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.instr.goes import _goes_get_chianti_em
    >>> longflux = u.Quantity([7e-6,7e-6], unit=u.W/u.m/u.m)
    >>> temp = u.Quantity([11, 11], unit=u.MK)
    >>> em = _goes_get_chianti_em(longflux, temp, satellite=15,
    ...                           abundances="coronal")  # doctest: +REMOTE_DATA
    >>> em  # doctest: +REMOTE_DATA
    <Quantity [3.45200672e+48, 3.45200672e+48] 1 / cm3>
    """
    # Check inputs are of correct type
    longflux = longflux.to(u.W/u.m**2)
    temp = temp.to(u.MK)
    # Ignore zero values raising a numpy warning here
    with np.errstate(invalid='ignore'):
        log10_temp = np.log10(temp.value)
    int(satellite)
    if satellite < 1:
        raise ValueError("satellite must be the number of a "
                         "valid GOES satellite (>1).")
    # if abundance input is valid create file suffix, abund, equalling
    # of 'cor' or 'pho'.
    if abundances == "coronal":
        data_file = manager.get('file_em_cor')
    elif abundances == "photospheric":
        data_file = manager.get('file_em_pho')
    else:
        raise ValueError("abundances must be a string equalling "
                         "'coronal' or 'photospheric'.")
    # check input arrays are of same length
    if len(longflux) != len(temp):
        raise ValueError("longflux and temp must have same number of "
                         "elements.")

    # Initialize lists to hold model data of temperature - long channel
    # flux relationship read in from csv file.
    modeltemp = []   # modelled temperature is in log_10 space in units of MK
    modelflux = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = f"longfluxGOES{satellite}"

    # Read data representing appropriate temperature--long flux
    # relationship depending on satellite number and assumed abundances.
    with open(data_file, "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelflux.append(float(row[label]))
    modeltemp = np.asarray(modeltemp)
    modelflux = np.asarray(modelflux)

    # Ensure input values of flux ratio are within limits of model table
    if np.min(log10_temp) < np.min(modeltemp) or \
       np.max(log10_temp) > np.max(modeltemp) or \
       np.isnan(np.min(log10_temp)):
        raise ValueError("All values in temp must be within the range "
                         "{} - {} MK.".format(np.min(10**modeltemp),
                                              np.max(10**modeltemp)))

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(log10_temp, spline, der=0)
    em = longflux.value/denom * 1e55
    em = u.Quantity(em, unit='cm**(-3)')

    return em


def calculate_radiative_loss_rate(goests, force_download=False,
                                  download_dir=None):
    """
    Calculates radiative loss rate from GOES observations.

    This function calculates the radiative loss rate as a function of
    time of solar soft X-ray-emitting plasma across all wavelengths given a
    LightCurve object containing GOES data.  The radiative loss rate is
    determined from the GOES isothermal temperature and volume emission
    measure as a function of time, as calculated by
    `~calculate_temperature_em()`.  See docstring of that function for more
    details.  If the LightCurve object does not contain the temperatures and
    emission measures, but only contain the GOES fluxes, then the temperature
    and emission measures are calculated using calculate_temperature_em().
    The unit of the resulting radiative loss rates is W.  Once
    the radiative loss rates have been found, they are returned as part of a
    new LightCurve object also containing the metadata, GOES fluxes and
    corresponding temperatures and emission measures of the input LightCurve
    object.

    Parameters
    ----------
    goests : `~sunpy.timeseries.sources.XRSTimeSeries`
        TimeSeries object containing GOES data.  The units of these
        data MUST be W/m^2 (flux), MK (temperature) and cm^-3
        (emission measure).  If LightCurve object does not contain
        temperature and emission measure values, they are calculated from
        the flux values using calculate_temperature_em().

    force_download : (optional) `bool`
        If True, the GOES radiative loss data file is downloaded even if
        already locally stored. It is important to do this if a new version
        of the file has been generated due to a new CHIANTI version being
        released or the launch of new GOES satellites.
        Default=False

    download_dir : (optional) `str`
        The directory to download the GOES radiative loss data file to.
        Default=SunPy default download directory

    Returns
    -------
    ts_new : `~sunpy.timeseries.sources.XRSTimeSeries`
        Contains same metadata and data as input LightCurve with the
        following additional data columns:

        | ts_new.to_dataframe().temperature - Array of temperature values [MK]
        | ts_new.to_dataframe().em - Array of volume emission measure values [cm**-3]
        | ts_new.to_dataframe().rad_loss_rate - radiative loss rate of the coronal soft
          X-ray-emitting plasma across all wavelengths [W]

    Notes
    -----
    The GOES radiative loss rates are calculated using a csv file containing
    a table of radiative loss rate per unit emission measure at various
    temperatures.  The appropriate values are then found via interpolation.
    This table was generated using CHIANTI atomic physics database employing
    the methods of [1]_.  Coronal abundances, a default
    density of 10**10 cm**-3, and ionization equilibrium of
    [2]_ were used.

    References
    ----------
    .. [1] Cox, D.P., Tucker, W.H. 1969, ApJ, 157, 1157, DOI: 10.1086/150144
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> import sunpy.timeseries as ts
    >>> from sunpy.instr.goes import calculate_radiative_loss_rate
    >>> from sunpy.data.sample import GOES_XRS_TIMESERIES  # doctest: +REMOTE_DATA
    >>> goests = ts.TimeSeries(GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> goests.to_dataframe()[0:10]  # doctest: +REMOTE_DATA
                                           xrsa          xrsb
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07
    >>> goests_new = calculate_radiative_loss_rate(goests)  # doctest: +REMOTE_DATA
    >>> goests_new.to_dataframe()[0:10]   # doctest:  +REMOTE_DATA
                                           xrsa          xrsb  temperature            em  rad_loss_rate
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07     3.503510  2.190626e+48   1.781001e+19
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07     3.534262  2.055847e+48   1.660031e+19
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07     3.518700  2.122771e+48   1.719931e+19
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07     3.550100  1.990333e+48   1.601718e+19
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07     3.518700  2.122771e+48   1.719931e+19
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07     3.550100  1.990333e+48   1.601718e+19
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07     3.518700  2.122771e+48   1.719931e+19
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07     3.534262  2.055847e+48   1.660031e+19
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07     3.550100  1.990333e+48   1.601718e+19
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07     3.534262  2.055847e+48   1.660031e+19
    """
    if not download_dir:
        download_dir = get_and_create_download_dir()
    # Check that input argument is of correct type
    if not isinstance(goests, timeseries.XRSTimeSeries):
        raise TypeError("goests must be a XRSTimeSeries object.")

    # extract temperature and emission measure from GOESLightCurve
    # object and change type to that required by _calc_rad_loss().
    # If LightCurve object does not contain temperature and
    # emission measure, calculate using calculate_temperature_em()
    if 'temperature' in goests.columns and 'em' in goests.columns:
        # Use copy.deepcopy for replicating meta and data so that input
        # lightcurve is not altered.
        ts_new = timeseries.XRSTimeSeries(meta=copy.deepcopy(goests.meta),
                                          data=copy.deepcopy(goests.to_dataframe()),
                                          units=copy.deepcopy(goests.units))
    else:
        ts_new = calculate_temperature_em(goests)
    temp = u.Quantity(np.asarray(ts_new.to_dataframe().temperature, dtype=np.float64),
                      unit=u.MK)
    em = u.Quantity(np.asarray(ts_new.to_dataframe().em, dtype=np.float64),
                    unit=u.cm**(-3))

    # Find radiative loss rate with _calc_rad_loss()
    rad_loss_out = _calc_rad_loss(temp, em, force_download=force_download,
                                  download_dir=download_dir)

    # Enter results into new version of GOES LightCurve Object
    ts_new = ts_new.add_column("rad_loss_rate", rad_loss_out['rad_loss_rate'].to("W"))

    return ts_new


@manager.require('file_rad_cor',
                 [urljoin(GOES_REMOTE_PATH, FILE_RAD_COR)],
                 'b56dccaa1035da46baa1a9251c4840107750d869de101d1811b506ceaec5828e')
@u.quantity_input
def _calc_rad_loss(temp: u.MK, em: u.cm**-3, obstime=None, force_download=False,
                   download_dir=None):
    """
    Finds radiative loss rate of coronal plasma over all wavelengths.

    This function calculates the radiative loss rate of solar coronal
    soft X-ray-emitting plasma across all wavelengths given an isothermal
    temperature and emission measure.  The units of the results are
    W.  This function is based on calc_rad_loss.pro in SSW IDL.
    In addition, if obstime keyword is set, giving the times to which
    the temperature and emission measure values correspond, the
    radiated losses integrated over time are also calculated.

    Parameters
    ----------
    temp : `~astropy.units.Quantity`
        Array containing the temperature of the coronal plasma at
        different times.  Units=[MK]

    em : `~astropy.units.Quantity`
        Array containing the emission measure of the coronal plasma
        at the same times corresponding to the temperatures in temp.
        Must be same length as temp.  Units=[cm**-3]

    obstime : (optional) array-like of `~sunpy.time.parse_time` parsable objects
        Array of measurement times to which temperature and
        emission measure values correspond.  Must be same length
        as temp and em.  If this keyword is set, the integrated
        radiated energy is calculated.

    force_download : (optional) bool
        If True, the GOES radiative loss data file is downloaded.  It is
        important to do this if a new version of the files has been
        generated due to a new CHIANTI version being released or the
        launch of new GOES satellites.
        Default=False

    download_dir : (optional) str
        The directory to download the GOES radiative loss data file to.
        Default=SunPy default download directory

    Returns
    -------
    rad_loss_out : `dict` of `~astropy.units.quantity.Quantity` objects
        Contains the following keys.

        | "rad_loss_rate" - radiative loss rate of the soft X-ray-emitting
           plasma across all wavelengths corresponding to temperatures and
           emission measures in temp and em Quantity inputs.
        | "rad_loss_cumul" - cumulative radiative losses as a function of
          time.  (Only if obstime kwarg is NOT None.)
        | "rad_loss_int" - total radiative losses as a function of time.
          (Only if obstime kwarg is not None.)  Array containing radiative
          loss rates of the coronal plasma corresponding to temperatures and
          emission measures in temp and em arrays.

    Notes
    -----
    This function calls a csv file containing a table of radiative loss
    rate per unit emission measure at various temperatures.  The
    appropriate values are then found via interpolation.  This table
    was generated using CHIANTI atomic physics database employing the
    methods of Cox & Tucker (1969).  Coronal abundances, a default
    density of 10**10 cm**-3, and ionization equilibrium of
    Mazzotta et al. (1998) were used.

    References
    ----------
    .. [1] Cox, D.P., Tucker, W.H. 1969, ApJ, 157, 1157, DOI: 10.1086/150144
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339, DOI: 10.1051/aas:1998330

    Examples
    --------
    >>> from sunpy.instr.goes import _calc_rad_loss
    >>> from astropy.units.quantity import Quantity
    >>> temp = Quantity([11.0, 11.0], unit="MK")
    >>> em = Quantity([4.0e+48, 4.0e+48], unit="cm**(-3)")
    >>> rad_loss = _calc_rad_loss(temp, em)  # doctest: +REMOTE_DATA
    >>> rad_loss["rad_loss_rate"]  # doctest: +REMOTE_DATA
    <Quantity [3.01851392e+19, 3.01851392e+19] J / s>
    """
    if not download_dir:
        download_dir = get_and_create_download_dir()
    # Check inputs are correct
    temp = temp.to(u.K)
    em = em.to(1/u.cm**3)
    if len(temp) != len(em):
        raise ValueError("temp and em must all have same number of elements.")

    # Initialize lists to hold model data of temperature - rad loss rate
    # relationship read in from csv file
    modeltemp = []    # modelled temperature is in log_10 space in units of MK
    model_loss_rate = []

    # Read data from csv file into lists, being sure to skip commented
    # lines beginning with "#"
    with open(manager.get('file_rad_cor'), "r") as csvfile:
        startline = csvfile.readlines()[7:]
        csvreader = csv.reader(startline, delimiter=" ")
        for row in csvreader:
            modeltemp.append(float(row[0]))
            model_loss_rate.append(float(row[1]))
    modeltemp = np.asarray(modeltemp)
    model_loss_rate = np.asarray(model_loss_rate)
    # Ensure input values of flux ratio are within limits of model table
    if temp.value.min() < modeltemp.min() or temp.value.max() > modeltemp.max():
        raise ValueError("All values in temp must be within the range " +
                         "{} - {} MK.".format(np.min(modeltemp/1e6),
                                              np.max(modeltemp/1e6)))
    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modeltemp, model_loss_rate, s=0)
    rad_loss = em.value * interpolate.splev(temp.value, spline, der=0)
    rad_loss = u.Quantity(rad_loss, unit='erg/s')
    rad_loss = rad_loss.to(u.J/u.s)

    # If obstime keyword giving measurement times is set, calculate
    # radiative losses integrated over time.
    if obstime is not None:
        # First ensure obstime is of same length as temp and em and of
        # correct type.
        n = len(temp)
        if len(obstime) != n:
            raise OSError("obstime must have same number of elements as "
                          "temp and em.")

        obstime = parse_time(obstime)

        # Check elements in obstime in chronological order
        _assert_chrono_order(obstime)
        # Next, get measurement times in seconds from time of first
        # measurement.
        obstime_seconds = (obstime - obstime[0]).sec
        # Finally, integrate using trapezoid rule
        rad_loss_int = trapz(rad_loss.value, obstime_seconds)
        rad_loss_int = u.Quantity(rad_loss_int, unit=rad_loss.unit*u.s)
        # Calculate cumulative radiated energy in each GOES channel as
        # a function of time.
        rad_loss_cumul = cumtrapz(rad_loss, obstime_seconds)
        rad_loss_cumul = u.Quantity(rad_loss_cumul, unit=rad_loss.unit*u.s)
        # Enter results into output dictionary.
        rad_loss_out = {"rad_loss_rate": rad_loss,
                        "rad_loss_cumul": rad_loss_cumul,
                        "rad_loss_int": rad_loss_int}
    else:
        rad_loss_out = {"rad_loss_rate": rad_loss}

    return rad_loss_out


def calculate_xray_luminosity(goests):
    """
    Calculates GOES solar X-ray luminosity.

    This function calculates the solar X-ray luminosity in the GOES
    wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based on the
    observed GOES fluxes.  The units of the results are W.  The calculation
    is made by simply assuming that the radiation is emitted isotropically,
    i.e. is distributed over a spherical surface area with a radius equal to
    the Sun-Earth distance.  Once the luminosity in each GOES passband is
    found, they are returned in a new LightCurve object also containing the
    metadata and data of the input LightCurve object.

    Parameters
    ----------
    goests : `~sunpy.timeseries.sources.XRSTimeSeries`
        LightCurve object containing GOES flux data which MUST
        be in units of W/m^2.

    Returns
    -------
    ts_new : `~sunpy.timeseries.sources.XRSTimeSeries`
        Contains same metadata and data as input LightCurve with the
        following additional data columns;

        | goests_new.to_dataframe().luminosity_xrsa - Xray luminosity in 0.5-4A channel
          unit=[W]
        | goests_new.to_dataframe().luminosity_xrsb - Xray luminosity in 1-8A channel
          unit=[W]

    Examples
    --------
    >>> import sunpy.timeseries as ts
    >>> from sunpy.instr.goes import calculate_xray_luminosity
    >>> from sunpy.data.sample import GOES_XRS_TIMESERIES  # doctest: +REMOTE_DATA
    >>> goests = ts.TimeSeries(GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> goests.to_dataframe()[0:10]  # doctest: +REMOTE_DATA
                                           xrsa          xrsb
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07
    >>> goests_new = calculate_xray_luminosity(goests)  # doctest: +REMOTE_DATA
    >>> goests_new.to_dataframe()[0:10]   # doctest:  +REMOTE_DATA
                                           xrsa          xrsb  luminosity_xrsa  luminosity_xrsb
    2011-06-06 23:59:59.961999893  1.000000e-09  1.887100e-07     2.896209e+14     5.465435e+16
    2011-06-07 00:00:02.008999944  1.000000e-09  1.834600e-07     2.896209e+14     5.313384e+16
    2011-06-07 00:00:04.058999896  1.000000e-09  1.860900e-07     2.896209e+14     5.389555e+16
    2011-06-07 00:00:06.104999900  1.000000e-09  1.808400e-07     2.896209e+14     5.237503e+16
    2011-06-07 00:00:08.151999950  1.000000e-09  1.860900e-07     2.896209e+14     5.389555e+16
    2011-06-07 00:00:10.201999903  1.000000e-09  1.808400e-07     2.896209e+14     5.237503e+16
    2011-06-07 00:00:12.248999953  1.000000e-09  1.860900e-07     2.896209e+14     5.389555e+16
    2011-06-07 00:00:14.298999906  1.000000e-09  1.834600e-07     2.896209e+14     5.313384e+16
    2011-06-07 00:00:16.344999909  1.000000e-09  1.808400e-07     2.896209e+14     5.237503e+16
    2011-06-07 00:00:18.391999960  1.000000e-09  1.834600e-07     2.896209e+14     5.313384e+16
    """
    # Check that input argument is of correct type
    if not isinstance(goests, timeseries.XRSTimeSeries):
        raise TypeError("goeslc must be a XRSTimeSeries object.")
    # Find temperature and emission measure with _goes_chianti_tem
    lx_out = _goes_lx(goests.quantity("xrsb"),
                      goests.quantity("xrsa"),
                      date=str(goests.to_dataframe().index[0]))
    # Enter results into new version of GOES LightCurve Object
    # Use copy.deepcopy for replicating meta and data so that input
    # lightcurve is not altered.
    ts_new = timeseries.XRSTimeSeries(meta=copy.deepcopy(goests.meta),
                                      data=copy.deepcopy(goests.to_dataframe()),
                                      units=copy.deepcopy(goests.units))
    ts_new = ts_new.add_column("luminosity_xrsa", lx_out["shortlum"].to("W"))
    ts_new = ts_new.add_column("luminosity_xrsb", lx_out["longlum"].to("W"))

    return ts_new


def _goes_lx(longflux, shortflux, obstime=None, date=None):
    """
    Calculates solar X-ray luminosity in GOES wavelength ranges.

    This function calculates the X-ray luminosity from the Sun in the
    GOES wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based
    on the observed GOES fluxes.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    longflux : `~astropy.units.Quantity`
        Array containing the observed GOES/XRS long channel flux.
        Units=[W/m**2]

    shortflux : `~astropy.units.Quantity`
        Array containing the observed GOES/XRS short channel flux.
        Units=[W/m**2]

    obstime : (optional) array-like of `~sunpy.time.parse_time` parsable objects
        Measurement times corresponding to each flux measurement.
        Assumes each pair of 0.5-4 and 1-8 angstrom flux measurements
        were taken simultaneously.

    date : (optional) `astropy.time.Time` object or valid date string.
        Date at which measurements were taken.  This is used to
        calculate the Sun-Earth distance.
        Default=None implies Sun-Earth distance is set to 1AU.

    Returns
    -------
    lx_out : `dict`
        dictionary containing the following fields.
        longlum : `~astropy.units.Quantity`
            Array of luminosity in the 1-8 angstroms range.

        shortlum : `~astropy.units.Quantity`
            Array of luminosity in the 0.5-4 angstroms range.

        longlum_int : (only present if obstime kwarg is set)

        shortlum_int : (only present if obstime kwarg is set)

    Notes
    -----
    This function calls _calc_xraylum() to calculate luminosities.
    For more information on how this is done, see docstring of that
    function.

    Examples
    --------
    >>> from sunpy.instr.goes import _goes_lx
    >>> from astropy.time import Time
    >>> from astropy.units.quantity import Quantity
    >>> longflux = Quantity([7e-6,7e-6,7e-6,7e-6,7e-6,7e-6], unit='W/m**2')
    >>> shortflux = Quantity([7e-7,7e-7,7e-7,7e-7,7e-7,7e-7], unit='W/m**2')
    >>> obstime = Time(['2014-1-1T0:0:0',
    ...                 '2014-1-1T0:0:2',
    ...                 '2014-1-1T0:0:4',
    ...                 '2014-1-1T0:0:6',
    ...                 '2014-1-1T0:0:8',
    ...                 '2014-1-1T0:0:10'])
    >>> lx_out = _goes_lx(longflux, shortflux, obstime)  # doctest: +REMOTE_DATA
    >>> lx_out["longlum"]  # doctest: +REMOTE_DATA
    <Quantity [1.96860565e+18, 1.96860565e+18, 1.96860565e+18, 1.96860565e+18,
                 1.96860565e+18, 1.96860565e+18] W>
    >>> lx_out["shortlum"]  # doctest: +REMOTE_DATA
    <Quantity [1.96860565e+17, 1.96860565e+17, 1.96860565e+17, 1.96860565e+17,
               1.96860565e+17, 1.96860565e+17] W>
    >>> lx_out["longlum_int"]  # doctest: +REMOTE_DATA
    <Quantity 1.96860565e+19 s W>
    >>> lx_out["shortlum_int"]  # doctest: +REMOTE_DATA
    <Quantity 1.96860565e+18 s W>
    """
    # Calculate X-ray luminosities
    longlum = _calc_xraylum(longflux, date=date)
    shortlum = _calc_xraylum(shortflux, date=date)

    # If obstime keyword giving measurement times is set, calculate
    # total energy radiated in the GOES bandpasses during the flare.
    if obstime is not None:
        # First ensure longflux, shortflux, and obstime are all of
        # equal length and obstime is of correct type.
        if not len(longflux) == len(shortflux) == len(obstime):
            raise ValueError("longflux, shortflux, and obstime must all have "
                             "same number of elements.")

        obstime = parse_time(obstime)

        # Check elements in obstime in chronological order
        _assert_chrono_order(obstime)

        # Next, get measurement times in seconds from time of first
        # measurement.
        obstime_seconds = (obstime - obstime[0]).sec

        # Finally, integrate using trapezoid rule
        longlum_int = trapz(longlum.value, obstime_seconds)
        longlum_int = u.Quantity(longlum_int, unit=longlum.unit*u.s)
        shortlum_int = trapz(shortlum.value, obstime_seconds)
        shortlum_int = u.Quantity(shortlum_int, unit=shortlum.unit*u.s)
        # Calculate cumulative radiated energy in each GOES channel as
        # a function of time.
        longlum_cumul = cumtrapz(longlum.value, obstime_seconds)
        longlum_cumul = u.Quantity(longlum_cumul, unit=longlum.unit*u.s)
        shortlum_cumul = cumtrapz(shortlum.value, obstime_seconds)
        shortlum_cumul = u.Quantity(shortlum_cumul,
                                    unit=shortlum.unit*u.s)
        lx_out = {"longlum": longlum, "shortlum": shortlum,
                  "longlum_cumul": longlum_cumul,
                  "shortlum_cumul": shortlum_cumul,
                  "longlum_int": longlum_int, "shortlum_int": shortlum_int}
    else:
        lx_out = {"longlum": longlum, "shortlum": shortlum}

    return lx_out


@u.quantity_input
def _calc_xraylum(flux: u.W/u.m/u.m, date=None):
    """
    Calculates solar luminosity based on observed flux observed at 1AU.

    This function calculates the luminosity from the Sun based
    on observed flux in W/m**2.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    flux : `~astropy.units.Quantity`
       Containing the observed solar flux.  Units=[W/m**2]

    date : (optional) `astropy.time.Time` object or valid date str
        Used to calculate a more accurate Sun-Earth distance based on
        Earth's orbit at that date.  If date is None, Sun-Earth
        distance is set to 1AU.

    Returns
    -------
    xraylum : `~astropy.units.Quantity` array with units=erg/s.
        Array of X-ray luminosity.

    Examples
    --------
    >>> from sunpy.instr.goes import _calc_xraylum
    >>> from astropy.units.quantity import Quantity
    >>> flux = Quantity([7e-6,7e-6], unit="W/m**2")
    >>> xraylum = _calc_xraylum(flux, date="2014-04-21")  # doctest: +REMOTE_DATA
    >>> xraylum  # doctest: +REMOTE_DATA
    <Quantity [1.98751663e+18, 1.98751663e+18] W>
    """
    if date is not None:
        date = parse_time(date)
        xraylum = 4 * np.pi * sun.earth_distance(date).to("m")**2 * flux
    else:
        xraylum = 4 * np.pi * constants.au.to("m")**2 * flux
    return xraylum


def flareclass_to_flux(flareclass):
    """
    Converts a GOES flare class into the corresponding X-ray flux.

    Parameters
    ----------
    flareclass : str
        The case-insensitive flare class (e.g., 'X3.2', 'm1.5', 'A9.6').

    Returns
    -------
    flux : `~astropy.units.Quantity`
        X-ray flux between 1 and 8 Angstroms as measured near Earth in W/m^2.

    Raises
    ------
    TypeError
        Input must be a string.

    Examples
    --------
    >>> from sunpy.instr.goes import flareclass_to_flux
    >>> flareclass_to_flux('A1.0')
    <Quantity 1.e-08 W / m2>
    >>> flareclass_to_flux('c4.7')
    <Quantity 4.7e-06 W / m2>
    >>> flareclass_to_flux('X2.4')
    <Quantity 0.00024 W / m2>
    """
    if not isinstance(flareclass, type('str')):
        raise TypeError("Input must be a string")
    # TODO should probably make sure the string is in the expected format.

    flareclass = flareclass.upper()
    # invert the conversion dictionary
    # conversion_dict = {v: k for k, v in GOES_CONVERSION_DICT.items()}
    return float(flareclass[1:]) * GOES_CONVERSION_DICT[flareclass[0]]


@u.quantity_input
def flux_to_flareclass(goesflux: u.watt/u.m**2):
    """
    Converts X-ray flux into the corresponding GOES flare class.

    Parameters
    ----------
    flux : `~astropy.units.Quantity`
        X-ray flux between 1 and 8 Angstroms (usually measured by GOES) as
        measured at the Earth in W/m^2

    Returns
    -------
    flareclass : str
        The flare class e.g.: 'X3.2', 'M1.5', 'A9.6'.

    Raises
    ------
    ValueError
        Flux cannot be negative.

    References
    ----------
    `Solar Flare Classification <https://en.wikipedia.org/wiki/Solar_flare#Classification>`_

    Examples
    --------
    >>> from sunpy.instr.goes import flux_to_flareclass
    >>> import astropy.units as u
    >>> flux_to_flareclass(1e-08 * u.watt/u.m**2)
    'A1'
    >>> flux_to_flareclass(4.7e-06 * u.watt/u.m**2)
    'C4.7'
    >>> flux_to_flareclass(0.00024 * u.watt/u.m**2)
    'X2.4'
    >>> flux_to_flareclass(7.8e-09 * u.watt/u.m**2)
    'A0.78'
    >>> flux_to_flareclass(0.00682 * u.watt/u.m**2)
    'X68.2'
    """

    if goesflux.value < 0:
        raise ValueError("Flux cannot be negative")

    decade = np.floor(np.log10(goesflux.to('W/m**2').value))
    # invert the conversion dictionary
    conversion_dict = {v: k for k, v in GOES_CONVERSION_DICT.items()}
    if decade < -8:
        str_class = "A"
        decade = -8
    elif decade > -4:
        str_class = "X"
        decade = -4
    else:
        str_class = conversion_dict.get(u.Quantity(10 ** decade, "W/m**2"))
    goes_subclass = 10 ** -decade * goesflux.to('W/m**2').value
    return f"{str_class}{goes_subclass:.3g}"


def _assert_chrono_order(obstime):
    chrono_check = obstime[1:] - obstime[:-1]
    if not all(val > TimeDelta(0*u.day) for val in chrono_check):
        raise ValueError("Elements of obstime must be in chronological order.")
