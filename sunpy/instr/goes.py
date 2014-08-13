from __future__ import absolute_import
from __future__ import division

import os.path
import datetime
import csv
import copy
import socket
from itertools import dropwhile

import numpy as np
from scipy import interpolate
from scipy.integrate import trapz
from scipy.integrate import cumtrapz
import astropy
from astropy import units
from astropy.units.quantity import Quantity
import pandas


from sunpy.net import hek
from sunpy.time import parse_time
from sunpy import config
from sunpy import lightcurve
from sunpy.util.net import check_download_file
from sunpy.sun import sun

__all__ = ['get_goes_event_list', 'temp_em', 'goes_chianti_tem']

# Check required data files are present in user's default download dir
# Define location where GOES data files are stored.
# Manually resolve the hostname
HOST = socket.gethostbyname_ex('hesperia.gsfc.nasa.gov')[-1][0]
GOES_REMOTE_PATH = "http://{0}/ssw/gen/idl/synoptic/goes/".format(HOST)
# Define location where data files should be downloaded to.
DATA_PATH = config.get("downloads", "download_dir")
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
    timerange: sunpy.time.TimeRange
        The time range to download the event list for.

    goes_class_filter: (optional) string
        A string specifying a minimum GOES class for inclusion in the list,
        e.g. 'M1', 'X2'.

    """
    # use HEK module to search for GOES events
    client = hek.HEKClient()
    event_type = 'FL'
    tstart = timerange.start()
    tend = timerange.end()

    # query the HEK for a list of events detected by the GOES instrument
    # between tstart and tend (using a GOES-class filter)
    if goes_class_filter:
        result = client.query(hek.attrs.Time(tstart, tend),
                              hek.attrs.EventType(event_type),
                              hek.attrs.FL.GOESCls > goes_class_filter,
                              hek.attrs.OBS.Observatory == 'GOES')
    else:
        result = client.query(hek.attrs.Time(tstart, tend),
                              hek.attrs.EventType(event_type),
                              hek.attrs.OBS.Observatory == 'GOES')

    # want to condense the results of the query into a more manageable
    # dictionary
    # keep event data, start time, peak time, end time, GOES-class,
    # location, active region source (as per GOES list standard)
    # make this into a list of dictionaries
    goes_event_list = []

    for r in result:
        goes_event = {
            'event_date': parse_time(r['event_starttime']).date().strftime(
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

def temp_em(goeslc, abundances="coronal",
            download=False, download_dir=DATA_PATH):
    """
    Calculates and adds temperature and EM to a GOESLightCurve.

    This function calculates the isothermal temperature and volume
    emission measure of the solar soft X-ray emitting plasma observed by
    the GOES/XRS.  This is done using the function goes_chianti_tem().
    See that function for more details.  Once the temperature and
    emission measure are found, they are added to a copy of the
    original GOESLightCurve object as goeslc.data.temperature and
    goeslc.data.em where goeslc is the GOESLightCurve object.

    Parameters
    ----------
    goeslc : GOESLightCurve object

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

    download_dir : (optional) string
        The directory to download the GOES temperature and emission measure
        data files to.
        Default=SunPy default download directory

    Returns
    -------
    goeslc.data.temperature : pandas.core.series.Series
        Array of temperature values [MK]

    goeslc.data.em : pandas.core.series.Series
        Array of volume emission measure values [10**49 cm**-3]

    Examples
    --------
    >>> import sunpy.lightcurve as lc
    >>> goeslc = lc.GOESLightCurve.create(time1, time2)
    >>> goeslc.data
                          xrsa   xrsb
    2014-01-01 00:00:00  7e-07  7e-06
    2014-01-01 00:00:02  7e-07  7e-06
    2014-01-01 00:00:04  7e-07  7e-06
    2014-01-01 00:00:06  7e-07  7e-06
    >>> goeslc_new = temp_em(goeslc)
    >>> goeslc_new.data
                          xrsa   xrsb  temperature              em
    2014-01-01 00:00:00  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:02  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:04  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:06  7e-07  7e-06  11.28295376  4.78577516e+48

    """
    # Check that input argument is of correct type
    if not isinstance(goeslc, lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")

    # Find temperature and emission measure with goes_chianti_tem
    temp, em = goes_chianti_tem(goeslc.data.xrsb, goeslc.data.xrsa,
                                satellite=goeslc.meta["TELESCOP"].split()[1],
                                date=goeslc.data.index[0],
                                abundances=abundances, download=download,
                                download_dir=download_dir)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["temperature"] = temp
    goeslc_new.data["em"] = em

    return goeslc_new

def goes_chianti_tem(longflux, shortflux, satellite=8,
                     date=datetime.datetime.today(), abundances="coronal",
                     download=False, download_dir=DATA_PATH):
    """
    Calculates temperature and emission measure from GOES/XRS data.

    This function calculates the isothermal temperature and volume
    emission measure of the solar soft X-ray emitting plasma observed by
    the GOES/XRS.  This is done using the observed flux ratio of the
    short (0.5-4 angstrom) to long (1-8 angstrom) channels.

    Parameters
    ----------
    longflux, shortflux :astropy Quantity
        Arrays containing the long and short GOES/XRS flux measurements
        respectively as a function of time.  Must be of same length. [W/m**2].

    satellite : int (optional)
        Number of GOES satellite used to make observations, important for
        correct calibration of data.
        Default=8

    date : datetime object or str
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

    download_dir : (optional) string
        The directory to download the GOES temperature and emission measure
        data files to.
        Default=SunPy default download directory

    Returns
    -------
    temp : astropy Quantity
        Array of temperature values of same length as longflux and
        shortflux. Units=[MK]

    em : astropy Quantity
        Array of volume emission measure values of same length as longflux
        and shortflux.  Units=[10**49 cm**-3]

    Notes
    -----
    The temperature and volume emission measure are calculated here
    using the methods of White et al. (2005) who used the
    CHIANTI atomic physics database to model the response of the ratio
    of the short (0.5-4 angstrom) to long (1-8 angstrom) channels of the
    XRSs onboard various GOES satellites.  This method assumes an
    isothermal plasma, the ionisation equilibria of
    Mazzotta et al. (1998), and a constant density of 10**10 cm**-3.
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
    Ratio values corresponding to suxh fluxes are set to 0.003.

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339

    Examples
    --------
    >>> longflux = np.array([7e-6, 7e-6])
    >>> shortflux = np.array([7e-7, 7e-7])
    >>> temp, em = goes_chianti_tem(longflux, shortflux, satellite=15,
                                    date='2014-04-16',
                                    abundances="coronal")
    >>> temp
    array([11.28295376, 11.28295376])
    >>> em
    array([  4.78577516e+48,   4.78577516e+48])

    """
    # ENSURE INPUTS ARE OF CORRECT TYPE AND VALID VALUES
    if type(longflux) != astropy.units.quantity.Quantity:
        longflux = Quantity(longflux, unit='W/m**2')
    if type(shortflux) != astropy.units.quantity.Quantity:
        shortflux = Quantity(shortflux, unit='W/m**2')
    int(satellite)
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
    if date < datetime.datetime(1983, 06, 28) and satellite == 6:
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
    index = np.logical_or(shortflux_corrected.value < 1e-10,
                          longflux_corrected.value < 3e-8)
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio.value[index] = 0.003

    # FIND TEMPERATURE AND EMISSION MEASURE FROM FUNCTIONS BELOW
    temp = _goes_get_chianti_temp(fluxratio, satellite=satellite,
                                  abundances=abundances, download=download,
                                  download_dir=download_dir)
    em = _goes_get_chianti_em(longflux_corrected, temp, satellite=satellite,
                              abundances=abundances, download=download,
                              download_dir=download_dir)
    return temp, em

def _goes_get_chianti_temp(fluxratio, satellite=8, abundances="coronal",
                           download=False, download_dir=DATA_PATH):
    """
    Calculates temperature from GOES flux ratio.

    This function calculates the isothermal temperature of the solar
    soft X-ray emitting plasma observed by the GOES/XRS from the
    observed flux ratio of the short (0.5-4 angstrom) to
    long (1-8 angstrom) channels.  This function is not intended to be
    called directly but by goes_chianti_tem(), although it can be used
    independently.  However, if used independently data preparation,
    such as correctly rescaling fluxes for some satellites etc. will
    not be carried out.  This is done in goes_chianti_tem().

    Parameters
    ----------
    fluxratio : astropy Quantity
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

    download_dir : (optional) string
        The directory to download the GOES temperature data file to.
        Default=SunPy default download directory

    Returns
    -------
    temp : numpy array
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

    These files were calculated using the methods of White et al. (2005)
    who used the CHIANTI atomic physics database to model the response
    of the ratio of the short (0.5-4 angstrom) to long (1-8 angstrom)
    channels of the XRSs onboard various GOES satellites.  This method
    assumes an isothermal plasma, the ionisation equilibria of
    Mazzotta et al. (1998), and a constant density of 10**10 cm**-3.
    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_get_chianti_temp.pro in
    SolarSoftWare written in IDL by Stephen White.

    For correct preparation of GOES data before calculating temperature
    see goes_chianti_tem() (Notes section of docstring).

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339

    Examples
    --------
    >>> fluxratio = np.array([0.1,0.1])
    >>> temp = _goes_get_chianti_temp(fluxratio, satellite=15,
                                    abundances="coronal")
    >>> temp
    array([11.28295376, 11.28295376])

    """
    # If download kwarg is True, or required data files cannot be
    # found locally, download required data files.
    check_download_file(FILE_TEMP_COR, GOES_REMOTE_PATH, download_dir,
                        replace=download)
    check_download_file(FILE_TEMP_PHO, GOES_REMOTE_PATH, download_dir,
                        replace=download)

    # check inputs are correct
    if type(fluxratio) != astropy.units.quantity.Quantity:
        fluxratio = Quantity(fluxratio)
    int(satellite)
    if satellite < 1:
        raise ValueError("satellite must be the number of a "
                         "valid GOES satellite (>1).")
    # if abundance input is valid create file suffix, abund, equalling
    # of 'cor' or 'pho'.
    if abundances == "coronal":
        data_file = FILE_TEMP_COR
    elif abundances == "photospheric":
        data_file = FILE_TEMP_PHO
    else:
        raise ValueError("abundances must be a string equalling "
                         "'coronal' or 'photospheric'.")

    # Initialize lists to hold model data of flux ratio - temperature
    # relationship read in from csv file
    modeltemp = [] # modelled temperature is in log_10 space in units of MK
    modelratio = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = "ratioGOES{0}".format(satellite)
    # Read data representing appropriate temperature--flux ratio
    # relationship depending on satellite number and assumed abundances.
    with open(os.path.join(DATA_PATH, data_file), "r") as csvfile:
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
    temp = Quantity(temp, unit='MK')

    return temp

def _goes_get_chianti_em(longflux, temp, satellite=8, abundances="coronal",
                         download=False, download_dir=DATA_PATH):
    """
    Calculates emission measure from GOES 1-8A flux and temperature.

    This function calculates the emission measure of the solar
    soft X-ray emitting plasma observed by the GOES/XRS from the
    the ratio of the isothermal temperature and observed long channel
    (1-8 angstrom) flux which scales with the emission measure.
    This function is not intended to be called directly but by
    goes_chianti_tem(), although it can be used independently.
    However, if used independently data preparation, such as correctly
    rescaling fluxes for some satellites etc. will not be carried out.
    This is done in goes_chianti_tem().

    Parameters
    ----------
    longflux : astropy Quantity
        Array containing the observed GOES/XRS long channel flux.
        Units=[W/m**2]

    temp : astropy Quantity
        Array containing the GOES temperature.  Units=[MK]

    satellite : int (optional)
        Number of GOES satellite used to make observations.
        Important for correct calibration of data.
        Default=8

    abundances : (optional) string equalling 'coronal' or 'photospheric'
        States whether photospheric or coronal abundances should be
        assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES emission measure data file is downloaded.
        It is important to do this if a new version of the file has been
        generated due to a new CHIANTI version being released or the launch of
        new GOES satellites since these file was last downloaded.
        Default=False

    download_dir : (optional) string
        The directory to download the GOES emission measure data file to.
        Default=SunPy default download directory

    Returns
    -------
    em : numpy array
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
    temperture of the emitting plasma for XRSs onboard various GOES
    satellites.  The emission measure can then be found by scaling the
    ratio of these two properties.  This method assumes an isothermal
    plasma, the ionisation equilibria of Mazzotta et al. (1998), and
    a constant density of 10**10 cm**-3.
    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_get_chianti_temp.pro in
    SolarSoftWare written in IDL by Stephen White.

    For correct preparation of GOES data before calculating temperature
    see goes_chianti_tem() (Notes section of docstring).

    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005,
        Sol. Phys., 227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., &
        Vittorio, N. 1998, A&AS, 133, 339

    Examples
    --------
    >>> longflux = np.array([7e-6,7e-6])
    >>> temp = np.array([11,11])
    >>> em = _goes_get_chianti_em(longflux, temp, satellite=15,
                                  abundances="coronal")
    >>> em
    array([  3.45200672e+48,   3.45200672e+48])

    """
    # If download kwarg is True, or required data files cannot be
    # found locally, download required data files.
    check_download_file(FILE_EM_COR, GOES_REMOTE_PATH, download_dir,
                        replace=download)
    check_download_file(FILE_EM_PHO, GOES_REMOTE_PATH, download_dir,
                        replace=download)

    # Check inputs are of correct type
    if type(longflux) != astropy.units.quantity.Quantity:
        longflux = Quantity(longflux, unit='W/m**2')
    else:
        longflux.to(units.W/units.m**2)
    if type(temp) != astropy.units.quantity.Quantity:
        temp = Quantity(temp, unit='MK')
    else:
        temp = temp.to(units.MK)
    int(satellite)
    if satellite < 1:
        raise ValueError("satellite must be the number of a "
                         "valid GOES satellite (>1).")
    # if abundance input is valid create file suffix, abund, equalling
    # of 'cor' or 'pho'.
    if abundances == "coronal":
        data_file = FILE_EM_COR
    elif abundances == "photospheric":
        data_file = FILE_EM_PHO
    else:
        raise ValueError("abundances must be a string equalling "
                         "'coronal' or 'photospheric'.")
    # check input arrays are of same length
    if len(longflux) != len(temp):
        raise ValueError("longflux and temp must have same number of "
                         "elements.")

    # Initialize lists to hold model data of temperature - long channel
    # flux relationship read in from csv file.
    modeltemp = [] # modelled temperature is in log_10 sapce in units of MK
    modelflux = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = "longfluxGOES{0}".format(satellite)

    # Read data representing appropriate temperature--long flux
    # relationship depending on satellite number and assumed abundances.
    with open(os.path.join(DATA_PATH, data_file), "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelflux.append(float(row[label]))
    modeltemp = np.asarray(modeltemp)
    modelflux = np.asarray(modelflux)

    # Ensure input values of flux ratio are within limits of model table
    if np.min(np.log10(temp.value)) < np.min(modeltemp) or \
      np.max(np.log10(temp.value)) > np.max(modeltemp) or \
      np.isnan(np.min(np.log10(temp.value))):
        raise ValueError("All values in temp must be within the range "
                         "{0} - {1} MK.".format(np.min(10**modeltemp),
                                                np.max(10**modeltemp)))

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(np.log10(temp.value), spline, der=0)
    em = longflux.value/denom * 1e55
    em = Quantity(em, unit='cm**(-3)')

    return em

def radiative_loss_rate(goeslc, force_download=False, download_dir=DATA_PATH):
    """
    Calculates flare radiative loss rate and adds it to GOESLightCurve.

    This function calculates the radiative loss rate as a function of
    time of solar coronal soft X-ray-emitting plasma across all
    wavelengths given a GOESLightCurve object.  The units of the
    results are erg/s. This is done by calling calc_rad_loss().
    For more information see documentation in that function.  Once
    the radiative loss rates have been found, it is added to a copy of
    the original GOESLightCurve object as goeslc_new.data.rad_loss_rate",
    where goeslc_new is the new GOESLightCurve object which is returned.

    Parameters
    ----------
    goeslc : GOESLightCurve object

    force_download : (optional) bool
        If True, the GOES radiative loss data file is downloaded.
        It is important to do this if a new version of the files has
        been generated due to a new CHIANTI version being released or
        the launch of new GOES satellites.
        Default=False

    download_dir : (optional) string
        The directory to download the GOES radiative loss data file to.
        Default=SunPy default download directory

    Returns
    -------
    goeslc_new : GOESLightCurve object
        a copy of the input GOESLightCurve object with an additional
        field, goeslc_new.data.rad_loss_rate
        (type=pandas.core.series.Series), which contains the radiative
        loss rate of the coronal soft X-ray-emitting plasma across all
        wavelengths in erg/s.  N.B. if the original GOESLightCurve
        object does not contain fields named goeslc_new.data.temperature
        and goeslc_new.data.em containing the temperature and emission
        measure values, these are also generated and added to goeslc_new
        using goes_chianti_tem() (See documentation for that function.)

    Examples
    --------
    >>> from sunpy.lightcurve as lc
    >>> goeslc = lc.GOESLightCurve.create(time1, time2)
    >>> goeslc.data
                               xrsa   xrsb
    2014-01-01 00:00:00  9.1873e-08  4e-06
    2014-01-01 00:00:02  9.1873e-08  4e-06
    2014-01-01 00:00:04  9.1873e-08  4e-06
    2014-01-01 00:00:06  9.2988e-08  4e-06
    2014-01-01 00:00:08  9.1873e-08  4e-06
    >>> goeslc_new = rad_loss_rate(goeslc)
    >>> goeslc_new.data
                        xrsa xrsb temperature            em  rad_loss_rate
    2014-01-01 00:00:00  ...  ...    6.270239  6.440648e+48 5.44914366e+26
    2014-01-01 00:00:02  ...  ...    6.270239  6.440648e+48 5.44914366e+26
    2014-01-01 00:00:04  ...  ...    6.273917  6.422207e+48 5.44914366e+26
    2014-01-01 00:00:06  ...  ...    6.304001  6.350370e+48 5.44914366e+26
    2014-01-01 00:00:08  ...  ...    6.277604  6.403795e+48 5.44914366e+26

    """
    # Check that input argument is of correct type
    if not isinstance(goeslc, lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")

    # extract temperature and emission measure from GOESLightCurve
    # object and change type to that required by calc_rad_loss().
    # If GOESLightCurve object does not contain temperature and
    # emission measure, calculate using temp_em()
    if 'temperature' in goeslc.data and 'em' in goeslc.data:
        goeslc_new = copy.deepcopy(goeslc)
    else:
        goeslc_new = temp_em(goeslc)
    temp = np.asarray(goeslc_new.data.temperature, dtype=np.float64)
    em = np.asarray(goeslc_new.data.em, dtype=np.float64)

    # Find radiative loss rate with calc_rad_loss()
    rad_loss_out = calc_rad_loss(temp, em, force_download=force_download,
                                 download_dir=download_dir)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new.data["rad_loss_rate"] = rad_loss_out["rad_loss_rate"]

    return goeslc_new

def calc_rad_loss(temp, em, obstime=None, cumulative=False,
                  force_download=False, download_dir=DATA_PATH):
    """
    Finds radiative loss rate of coronal plasma over all wavelengths.

    This function calculates the radiative loss rate of solar coronal
    soft X-ray-emitting plasma across all wavelengths given an isothermal
    temperature and emission measure.  The units of the results are
    erg/s.  This function is based on calc_rad_loss.pro in SSW IDL.
    In addition, if obstime keyword is set, giving the times to which
    the temperature and emission measure values correspond, the
    radiated losses integrated over time are also calculated.

    Parameters
    ----------
    temp : astropy Quantity
        Array containing the temperature of the coronal plasma at
        different times.  Units=[MK]

    em : astropy Quantity
        Array containing the emission measure of the coronal plasma
        at the same times corresponding to the temperatures in temp.
        Must be same length as temp.  Units=[cm**-3]

    obstime : (optional) array-like of datetime objects
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

    download_dir : (optional) string
        The directory to download the GOES radiative loss data file to.
        Default=SunPy default download directory

    Returns
    -------
    radlossrate : numpy ndarray, dtype=float, units=[erg]
        Array containing radiative loss rates of the coronal plasma
        corresponding to temperatures and emission measures in temp and
        em arrays.

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
    .. [1] Cox, D. P., Tucker, W. H. 1969, ApJ, 157, 1157
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339

    Examples
    --------
    >>> temp = np.array([11.0, 11.0])
    >>> em = np.array([4.0e+48, 4.0e+48])
    >>> rad_loss = calc_rad_loss(temp, em)
    >>> rad_loss["rad_loss_rate"]
    array([  3.01851392e+26,   3.01851392e+26])
    """
    # Check inputs are correct
    if type(temp) != astropy.units.quantity.Quantity:
        temp = Quantity(temp, unit='MK')
    temp = temp.to(units.K)
    if type(em) != astropy.units.quantity.Quantity:
        em = Quantity(em, unit='cm**(-3)')
    else:
        em = em.to(1/units.cm**3)
    if len(temp) != len(em):
        raise ValueError("temp and em must all have same number of elements.")
    # If force_download kwarg is True, or required data files cannot be
    # found locally, download required data files.
    check_download_file(FILE_RAD_COR, GOES_REMOTE_PATH, download_dir,
                        replace=force_download)

    # Initialize lists to hold model data of temperature - rad loss rate
    # relationship read in from csv file
    modeltemp = [] # modelled temperature is in log_10 sapce in units of MK
    model_loss_rate = []

    # Read data from csv file into lists, being sure to skip commented
    # lines begining with "#"
    with open(os.path.join(DATA_PATH, FILE_RAD_COR),
              "r") as csvfile:
        startline = csvfile.readlines()[7:]
        csvreader = csv.reader(startline, delimiter=" ")
        for row in csvreader:
            modeltemp.append(float(row[0]))
            model_loss_rate.append(float(row[1]))
    modeltemp = np.asarray(modeltemp)
    model_loss_rate = np.asarray(model_loss_rate)
    # Ensure input values of flux ratio are within limits of model table
    if np.min(temp.value) < np.min(modeltemp) or \
        np.max(temp.value) > np.max(modeltemp):
        raise ValueError("All values in temp must be within the range " +
                         "{0} - {1} MK.".format(np.min(modeltemp/1e6),
                                                np.max(modeltemp/1e6)))
    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modeltemp, model_loss_rate, s=0)
    rad_loss = em.value * interpolate.splev(temp.value, spline, der=0)
    rad_loss = Quantity(rad_loss, unit='erg/s')
    rad_loss = rad_loss.to(units.J/units.s)

    # If obstime keyword giving measurement times is set, calculate
    # radiative losses intergrated over time.
    if obstime is not None:
        # First ensure obstime is of same length as temp and em and of
        # correct type.
        n = len(temp)
        if len(obstime) != n:
            raise IOError("obstime must have same number of elements as "
                          "temp and em.")
        if type(obstime) == pandas.tseries.index.DatetimeIndex:
            obstime = obstime.to_pydatetime
        if any(type(obst) == str for obst in obstime):
            parse_time(obstime)
        if not all(type(obst) == datetime.datetime for obst in obstime):
            raise TypeError("obstime must be an array-like whose elements are"
                            " convertible to datetime objects.")
        # Next, get measurement times in seconds from time of first
        # measurement.
        obstime_seconds = np.array([(ot-obstime[0]).total_seconds()
                                    for ot in obstime], dtype="float64")
        # Finally, integrate using trapezoid rule
        rad_loss_int = trapz(rad_loss.value, obstime_seconds)
        rad_loss_int = Quantity(rad_loss_int, unit='J')
        # If cumulative kwarg True, calculate cumulative radiated energy
        # in each GOES channel as a function of time.
        if cumulative:
            rad_loss_cumul = cumtrapz(rad_loss, obstime_seconds)
            rad_loss_cumul = Quantity(rad_loss_cumul, unit='J')
            # Enter results into output dictionary.
            rad_loss_out = {"rad_loss_rate":rad_loss,
                            "rad_loss_cumul" : rad_loss_cumul,
                            "rad_loss_int":rad_loss_int}
        else:
            rad_loss_out = {"rad_loss_rate":rad_loss,
                            "rad_loss_int":rad_loss_int}
    else:
        # Ensure cumulative kwarg wasn't set without setting obstime.
        if cumulative:
            raise ValueError("cumulative keyword is True but obstime keyword "
                             "is None.  In order to calculate cumulative "
                             "radiated losses, cumulative must be True and "
                             "measurement times must be given via the "
                             "obstime kwarg.")
        # If keyword assignments are OK, enter results into output
        # dictionary.
        rad_loss_out = {"rad_loss_rate":rad_loss}

    return rad_loss_out

def xray_luminosity(goeslc):
    """
    Calculates and adds GOES solar X-ray luminosity to a GOESLightCurve.

    This function calculates the solar X-ray luminosity in the GOES
    wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based on the
    observed GOES fluxes.  The units of the results are erg/s. This is
    done by calling goes_lx().  This function assumes that the
    radiation is emitted isotropically, i.e. is distributed over a
    spherical surface area with a radius equal to the Sun-Earth
    distance.  Once the luminosity in each GOES passband is found, they
    are added to a copy of the original GOESLightCurve object as
    goeslc_new.data.luminosity_xrsa (for the 0.5-4 angstrom channel) and
    goeslc_new.data.luminosity_xrsb (for the 1-8 angstrom channel),
    where goeslc_new is the new GOESLightCurve object which is returned.

    Parameters
    ----------
    goeslc : GOESLightCurve object

    Returns
    -------
    goeslc_new : GOESLightCurve object
        A copy of the input GOESLightCurve object with two additional
        fields, goeslc_new.data.luminosity_xrsa and
        goeslc_new.data.luminosity_xrsb (each of type
        pandas.core.series.Series) which hold the X-ray luminosity in
        the 0.5-4 and 1-8 angstrom wavelength ranges, respectively.
        Units=[erg/s].

    Examples
    --------
    >>> from sunpy.lightcurve as lc
    >>> goeslc = lc.GOESLightCurve.create(time1, time2)
    >>> goeslc.data
                          xrsa   xrsb
    2014-01-01 00:00:00  7e-07  7e-06
    2014-01-01 00:00:02  7e-07  7e-06
    2014-01-01 00:00:04  7e-07  7e-06
    2014-01-01 00:00:06  7e-07  7e-06
    >>> goeslc_new = xray_luminosity(goeslc)
    >>> goeslc_new.data
                          xrsa   xrsb    luminosity_xrsa luminosity_xrsb
    2014-01-01 00:00:00  7e-07  7e-06     1.96860565e+24  1.96860565e+25
    2014-01-01 00:00:02  7e-07  7e-06     1.96860565e+24  1.96860565e+25
    2014-01-01 00:00:04  7e-07  7e-06     1.96860565e+24  1.96860565e+25
    2014-01-01 00:00:06  7e-07  7e-06     1.96860565e+24  1.96860565e+25

    """
    # Check that input argument is of correct type
    if not isinstance(goeslc, lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")
    # Find temperature and emission measure with goes_chianti_tem
    lx_out = goes_lx(goeslc.data.xrsb, goeslc.data.xrsa,
                     date=str(goeslc.data.index[0]))
    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["luminosity_xrsa"] = lx_out["shortlum"]
    goeslc_new.data["luminosity_xrsb"] = lx_out["longlum"]

    return goeslc_new

def goes_lx(longflux, shortflux, obstime=None, date=None, cumulative=False):
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
    longflux : astropy Quantity
        Array containing the observed GOES/XRS long channel flux.
        Units=[W/m**2]

    shortflux : astropy Quantity
        Array containing the observed GOES/XRS short channel flux.
        Units=[W/m**2]

    obstime : (optional) array-like of datetime objects
        Measurement times corresponding to each flux measurement.
        Assumes each pair of 0.5-4 and 1-8 angstrom flux measurements
        were taken simultaneously.

    date : (optional) datetime object or valid date string.
        Date at which measurements were taken.

    cumulative : (optional) bool
        If True and obstime is set, the cumulative radiated energy in
        each of the GOES wavelength bands is calculated as a function
        of time.
        Default=False

    Returns
    -------
    lx_out : dictionary
        dictionary containing the following fields.
        longlum : astropy Quantity
            Array of luminosity in the 1-8 angstroms range.

        shortlum : astropy Quantity
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
    >>> from datetime import datetime
    >>> from astropy.units.quantity import Quantity
    >>> longflux = Quantity([7e-6,7e-6,7e-6,7e-6,7e-6,7e-6], unit='W/m**2')
    >>> shortflux = Quantity([7e-7,7e-7,7e-7,7e-7,7e-7,7e-7], unit='W/m**2')
    >>> obstime = np.array([datetime(2014,1,1,0,0,0),
                            datetime(2014,1,1,0,0,2),
                            datetime(2014,1,1,0,0,4),
                            datetime(2014,1,1,0,0,6),
                            datetime(2014,1,1,0,0,8),
                            datetime(2014,1,1,0,0,10),], dtype=object)
    >>> lx_out = goes_lx(longflux, shortflux, obstime)
    >>> lx_out["longlum"]
    array([  1.96860565e+25,   1.96860565e+25,   1.96860565e+25,
             1.96860565e+25,   1.96860565e+25,   1.96860565e+25])
    >>> lx_out["shortlum"]
    array([  1.96860565e+24,   1.96860565e+24,   1.96860565e+24,
             1.96860565e+24,   1.96860565e+24,   1.96860565e+24])
    >>> lx_out["longlum_int"]
    1.96860565412e+26
    >>> lx_out["shortlum_int"]
    1.96860565412e+25

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
        if type(obstime) == pandas.tseries.index.DatetimeIndex:
            obstime = obstime.to_pydatetime
        if any(type(obst) == str for obst in obstime):
            parse_time(obstime)
        if not all(type(obst) == datetime.datetime for obst in obstime):
            raise TypeError("obstime must be an array-like whose elements are"
                            " convertible to datetime objects.")
        # Next, get measurement times in seconds from time of first
        # measurement.
        obstime_seconds = np.array([(ot-obstime[0]).total_seconds()
                                    for ot in obstime], dtype="float64")
        # Finally, integrate using trapezoid rule
        longlum_int = trapz(longlum.value, obstime_seconds)
        longlum_int = Quantity(longlum_int, unit=longlum.unit*units.s)
        longlum_int = longlum_int.to(units.J)
        shortlum_int = trapz(shortlum.value, obstime_seconds)
        shortlum_int = Quantity(shortlum_int, unit=shortlum.unit*units.s)
        shortlum_int = shortlum_int.to(units.J)
        # If cumulative kwarg True, calculate cumulative radiated energy
        # in each GOES channel as a function of time.
        if cumulative is True:
            longlum_cumul = cumtrapz(longlum.value, obstime_seconds)
            longlum_cumul = Quantity(longlum_cumul, unit='J/s')
            shortlum_cumul = cumtrapz(shortlum.value, obstime_seconds)
            shortlum_cumul = Quantity(shortlum_cumul, unit='J/s')
            lx_out = {"longlum":longlum, "shortlum":shortlum,
                      "longlum_cumul":longlum_cumul,
                      "shortlum_cumul":shortlum_cumul,
                      "longlum_int":longlum_int, "shortlum_int":shortlum_int}
        else:
            lx_out = {"longlum":longlum, "shortlum":shortlum,
                      "longlum_int":longlum_int, "shortlum_int":shortlum_int}
    else:
        # Ensure cumulative kwarg wasn't set without setting obstime.
        if cumulative is True:
            raise IOError("cumulative keyword is True but obstime keyword is "
                          "None.  In order to calculate cumulative X-ray "
                          "radiated energies, cumulative must be True and "
                          "measurement times must be given via the obstime "
                          "keyword.")
        # If keyword assignments are OK, enter results into output
        # dictionary.
        lx_out = {"longlum":longlum, "shortlum":shortlum}

    return lx_out

def _calc_xraylum(flux, date=None):
    """
    Calculates solar luminosity based on observed flux observed at 1AU.

    This function calculates the luminosity from the Sun based
    on observed flux in W/m**2.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    flux : astropy Quantity
       Containing the observed solar flux.  Units=[W/m**2]

    date : (optional) datetime object or valid date string
        Used to calculate a more accurate Sun-Earth distance based on
        Earth's orbit at that date.  If date is not set, standard value
        for 1AU used.

    Returns
    -------
    xraylum : numpy array, dtype=float, units=erg/s.
        Array of X-ray luminosity.

    Notes
    -----
    To convert from W/m**2 to erg/s:
    1 W = 1 J/s = 10**7 erg/s
    1 W/m**2 = 4*pi * AU**2 * 10**7 erg/s, where AU is the Sun-Earth
    distance in metres.

    Examples
    --------
    >>> flux = np.array([7e-6,7e-6])
    >>> xraylum = _calc_xraylum(flux, date="2014-04-21")
    >>> xraylum
    array([  1.98649103e+25,   1.98649103e+25])

    """
    # Ensure input is of correct type
    if type(flux) != astropy.units.quantity.Quantity:
        flux = Quantity(flux, unit='W/m**2')
    else:
        flux = flux.to(units.W/units.m**2)
    if date is not None:
        date = parse_time(date)
        xraylum = 4 * np.pi * (sun.sunearth_distance(t=date))**2 * flux
    else:
        xraylum = 4 * np.pi * (sun.sunearth_distance())**2 * flux
    xraylum = xraylum.si
    return xraylum
