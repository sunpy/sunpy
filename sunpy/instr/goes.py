from __future__ import absolute_import
from __future__ import division

import os.path
import datetime
import csv
import copy
import urllib
import urllib2
from itertools import dropwhile

import numpy as np
import scipy.interpolate as interpolate

from sunpy.net import hek
from sunpy.time import parse_time
from sunpy.sun import sun
from sunpy import config
import sunpy.lightcurve

__all__ = ['get_goes_event_list', 'temp_em', 'goes_chianti_tem']

# Check required data files are present in user's default download dir
# Define location where GOES data files are stored.
GOES_REMOTE_PATH = "http://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/"
# Define location where data files should be downloaded to.
DATA_PATH = config.get("downloads", "download_dir")
# Define variables for file names
FILE_TEMP_COR = "goes_chianti_temp_cor.csv"
FILE_TEMP_PHO = "goes_chianti_temp_pho.csv"
FILE_EM_COR = "goes_chianti_em_cor.csv"
FILE_EM_PHO = "goes_chianti_em_pho.csv"

def get_goes_event_list(trange, goes_class_filter=None):
    """
    Retrieve list of flares detected by GOES within a given time range.

    Parameters
    ----------
    trange: a SunPy TimeRange object

    goes_class_filter: (optional) string
                       a string specifying a minimum GOES class for
                       inclusion in the list, e.g. 'M1', 'X2'.

    """
    # use HEK module to search for GOES events
    client = hek.HEKClient()
    event_type = 'FL'
    tstart = trange.start()
    tend = trange.end()

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
            'event_date': parse_time(r['event_starttime']).date().strftime('%Y-%m-%d'),
            'start_time': parse_time(r['event_starttime']),
            'peak_time': parse_time(r['event_peaktime']),
            'end_time': parse_time(r['event_endtime']),
            'goes_class': str(r['fl_goescls']),
            'goes_location': (r['event_coord1'], r['event_coord2']),
            'noaa_active_region': r['ar_noaanum']
            }
        goes_event_list.append(goes_event)

    return goes_event_list

def temp_em(goeslc, abundances="coronal", download=False):
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
    abundances : (optional) string equalling either 'coronal' or
                 'photospheric'.
                 States whether photospheric or coronal abundances
                 should be assumed.
                 Default='coronal'
    download : (optional) bool
               If True, the GOES temperature and emission measure data
               files are downloaded.  It is important to do this if a
               new version of the files has been generated due to a new
               CHIANTI version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

    Returns
    -------
    goeslc.data.temperature : pandas.core.series.Series
                              Array of temperature values [MK]
    goeslc.data.em : pandas.core.series.Series
                     Array of volume emission measure values
                     [10**49 cm**-3]

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
    >>> goeslc_new = temp_em(goeslc)
    >>> goeslc_new.data
                          xrsa   xrsb  temperature              em
    2014-01-01 00:00:00  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:02  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:04  7e-07  7e-06  11.28295376  4.78577516e+48
    2014-01-01 00:00:06  7e-07  7e-06  11.28295376  4.78577516e+48

    """

    # Check that input argument is of correct type
    if not isinstance(goeslc, sunpy.lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")

    # Find temperature and emission measure with goes_chianti_tem
    temp, em = goes_chianti_tem(goeslc.data.xrsb, goeslc.data.xrsa,
                                satellite=goeslc.meta["TELESCOP"].split()[1],
                                date=goeslc.data.index[0],
                                abundances=abundances, download=download)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["temperature"] = temp
    goeslc_new.data["em"] = em

    return goeslc_new

def goes_chianti_tem(longflux, shortflux, satellite=8,
                     date=datetime.datetime.today(), abundances="coronal",
                     download=False):
    """
    Calculates temperature and emission measure from GOES/XRS data.

    This function calculates the isothermal temperature and volume
    emission measure of the solar soft X-ray emitting plasma observed by
    the GOES/XRS.  This is done using the observed flux ratio of the
    short (0.5-4 angstrom) to long (1-8 angstrom) channels.

    Parameters
    ----------
    longflux, shortflux : ndarray or array-like which can be converted
                          to float64 type, such as an np.array, tuple,
                          list.
                          Arrays containing the long and short GOES/XRS
                          flux measurements respectively as a function
                          of time.  Must be of same length. [W/m**2].
    satellite : int (optional)
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    date : datetime object or str
           Date when observations made.  Important for correct
           calibration.  Default=today
    abundances : (optional) string equalling either 'coronal' or
                 'photospheric'.
                 States whether photospheric or coronal abundances
                 should be assumed.
                 Default='coronal'
    download : (optional) bool
               If True, the GOES temperature and emission measure data
               files are downloaded.  It is important to do this if a
               new version of the files has been generated due to a new
               CHIANTI version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

    Returns
    -------
    temp : numpy array
           Array of temperature values of same length as longflux and
           shortflux.  [MK]
    em : numpy array
         Array of volume emission measure values of same length as
         longflux and shortflux.  [10**49 cm**-3]

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
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005, Sol. Phys.,
       227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339

    Examples
    --------
    >>> longflux = np.array([7e-6, 7e-6])
    >>> shortflux = np.array([7e-7, 7e-7])
    >>> temp, em = goes_chianti_tem(longflux, shortflux, satellite=15,
                                    date='2014-04-16', abundances="coronal")
    >>> temp
    array([11.28295376, 11.28295376])
    >>> em
    array([  4.78577516e+48,   4.78577516e+48])

    """
    # ENSURE INPUTS ARE OF CORRECT TYPE AND VALID VALUES
    longflux = np.asanyarray(longflux, dtype=np.float64)
    shortflux = np.asanyarray(shortflux, dtype=np.float64)
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
        longflux_corrected = longflux * (4.43/5.32)
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
    index = np.logical_or(shortflux_corrected < 1e-10,
                          longflux_corrected < 3e-8)
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio[index] = 0.003

    # FIND TEMPERATURE AND EMISSION MEASURE FROM FUNCTIONS BELOW
    temp = _goes_get_chianti_temp(fluxratio, satellite=satellite,
                                  abundances=abundances, download=download)
    em = _goes_get_chianti_em(longflux_corrected, temp, satellite=satellite,
                              abundances=abundances, download=download)
    return temp, em

def _goes_get_chianti_temp(fluxratio, satellite=8, abundances="coronal",
                           download=False):
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
    fluxratio : ndarray or array-like which can be converted to float64
                type, such as an np.array, tuple, list.
                Array containing the ratio of short channel to long
                channel GOES/XRS flux measurements.
    satellite : int (optional)
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    abundances : (optional) string equalling either 'coronal' or
                 'photospheric'.
                 States whether photospheric or coronal abundances
                 should be assumed.
                 Default='coronal'
    download : (optional) bool
               If True, the GOES temperature data files are
               downloaded.  It is important to do this if a new version
               of the files has been generated due to a new CHIANTI
               version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

    Returns
    -------
    temp : numpy array
           Array of temperature values of same length as longflux and
           shortflux.  [MK]

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
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005, Sol. Phys.,
       227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339

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
    _check_download_file(FILE_TEMP_COR, GOES_REMOTE_PATH, localpath=DATA_PATH,
                         force_download=download)
    _check_download_file(FILE_TEMP_PHO, GOES_REMOTE_PATH, localpath=DATA_PATH,
                         force_download=download)

    # check inputs are correct
    fluxratio = np.asanyarray(fluxratio, dtype=np.float64)
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
    temp = 10.**interpolate.splev(fluxratio, spline, der=0)

    return temp

def _goes_get_chianti_em(longflux, temp, satellite=8, abundances="coronal",
                         download=False):
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
    longflux : ndarray or array-like which can be converted to float64
               type, such as an np.array, tuple, list.
               Array containing the observed GOES/XRS long channel flux
    temp : ndarray or array-like which can be converted to float64
           type, such as an np.array, tuple, list.
           Array containing the GOES temperature
    satellite : int (optional)
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    abundances : (optional) string equalling either 'coronal' or
                 'photospheric'.
                 States whether photospheric or coronal abundances
                 should be assumed.
                 Default='coronal'
    download : (optional) bool
               If True, the GOES emission measure data files are
               downloaded.  It is important to do this if a new version
               of the files has been generated due to a new CHIANTI
               version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

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
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005, Sol. Phys.,
       227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339

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
    _check_download_file(FILE_EM_COR, GOES_REMOTE_PATH, localpath=DATA_PATH,
                         force_download=download)
    _check_download_file(FILE_EM_PHO, GOES_REMOTE_PATH, localpath=DATA_PATH,
                         force_download=download)

    # If download kwarg is True, download required data files
    if download:
        urllib.urlretrieve(os.path.join(GOES_REMOTE_PATH, FILE_EM_COR),
                           os.path.join(localpath, FILE_EM_COR))
        urllib.urlretrieve(os.path.join(GOES_REMOTE_PATH, FILE_EM_PHO),
                           os.path.join(localpath, FILE_EM_PHO))

    # Check inputs are of correct type
    longflux = np.asanyarray(longflux, dtype=np.float64)
    temp = np.asanyarray(temp, dtype=np.float64)
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
    if np.min(np.log10(temp)) < np.min(modeltemp) or \
      np.max(np.log10(temp)) > np.max(modeltemp):
        raise ValueError("All values in temp must be within the range "
                         "{0} - {1} MK.".format(np.min(10**modeltemp),
                                                np.max(10**modeltemp)))

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(np.log10(temp), spline, der=0)
    em = longflux/denom * 1e55

    return em

def rad_loss_rate(goeslc, download=False):
    """
    Calculates and adds solar radiative loss rate to a GOESLightCurve.

    This function calculates the radiative loss rate as a function of
    time of solar coronal soft X-ray-emitting plasma across all
    wavelengths given a GOESLightCurve object.  The units of the
    results are erg/s. This is done by calling calc_rad_loss().
    For more information see documentation in that function.  Once
    the radiative loss rates have been found, it is added to a copy of
    the original GOESLightCurve object as goeslc.data.rad_loss_rate",
    where goeslc is the GOESLightCurve object.

    Parameters
    ----------
    goeslc : GOESLightCurve object
    download : (optional) bool
               If True, the GOES CHIANTI radiative loss data
               files are downloaded.  It is important to do this if a
               new version of the files has been generated due to a new
               CHIANTI version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

    Returns
    -------
    goeslc.data.rad_loss_rate : pandas.core.series.Series
                                Array of radiative loss rate of the
                                coronal soft X-ray-emitting plasma
                                across all wavelengths. [erg/s]

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
    >>> goeslc_new = rad_loss_rate(goeslc)
    >>> goeslc_new.data
                          xrsa   xrsb  luminosity_xrsa  luminosity_xrsb
    2014-01-01 00:00:00  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:02  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:04  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:06  7e-07  7e-06     1.903523e+24     1.903523e+25

    """

    # Check that input argument is of correct type
    if not isinstance(goeslc, sunpy.lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")

    # extract temperature and emission measure from GOESLightCurve
    # object and change type to that required by calc_rad_loss().
    # If GOESLightCurve object does not contain temperature and
    # emission measure, calculate using temp_em()
    try:
        temp = np.asarray(goeslc.data.temperature, dtype=np.float64)
        em = np.asarray(goeslc.data.em, dtype=np.float64)
    except AttributeError as error:
        if error.message ==
            "'DataFrame' object has no attribute 'temperature'" or
            error.message == "'DataFrame' object has no attribute 'em'":
            goeslc_new = temp_em(goeslc)
            temp = np.asarray(goeslc_new.data.temperature, dtype=np.float64)
            em = np.asarray(goeslc_new.data.em, dtype=np.float64)
        else:
            raise error
    else:
        goeslc_new = copy.deepcopy(goeslc)

    # Find radiative loss rate with calc_rad_loss()
    rad_loss_out = calc_rad_loss(temp, em)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new.data["rad_loss_rate"] = rad_loss_out["rad_loss_rate"]

    return goeslc_new

def calc_rad_loss(temp, em, obstime=None, Download=False):
    """
    Finds radiative loss rate of solar SXR plasma over all wavelengths.

    This function calculates the luminosity of solar coronal soft
    X-ray-emitting plasma across all wavelengths given an isothermal
    temperature and emission measure.  The units of the results are
    erg/s.  This function is based on calc_rad_loss.pro in SSW IDL.
    In addition, if obstime keyword is set, giving the times to which
    the temperature and emission measure values correspond, the
    radiated losses integrated over time are also calculated.

    Parameters
    ----------
    temp : ndarray or array-like which can be converted to float64
           type, such as an np.array, tuple, list.  Units=[MK]
           Array containing the temperature of the coronal plasma at
           different times.
    em : ndarray or array-like which can be converted to float64 type,
         such as an np.array, tuple, list.  Units=[cm**-3]
         Array containing the emission measure of the coronal plasma
         at the same times corresponding to the temperatures in temp.
         Must be same length as temp
    obstime : (optional) ndaray array or array-like whose entries are
              (or can be converted to) datetime64 type, e.g. np.array,
              list, string array.
              Array of measurement times to which temperature and
              emission measure values correspond.  Must be same length
              as temp and em.  If this keyword is set, the integrated
              radiated energy is calculated.
    download : (optional) bool
               If True, the GOES temperature data files are
               downloaded.  It is important to do this if a new version
               of the files has been generated due to a new CHIANTI
               version being released or the launch of new GOES
               satellites since these files were originally downloaded.
               Default=False

    Returns
    -------
    radlossrate : numpy ndarray, dtype=float, units=[erg]
                  Array containing radiative loss rates of the coronal
                  plasma corresponding to temperatures and emission
                  measures in temp and em arrays.

    Notes
    -----
    This function calls a csv file containing a table of radiative loss
    rate per unit emission measure at various temperatures.  The
    appropriate values are then found via interpolation.
    This table was generated using CHIANTI atomic physics database
    employing the methods of Cox & Tucker (1969).  Coronal abundances,
    default density of 10**10 cm**-3, and ionization equilibrium of
    Mazzotta et al. (1998) were used.

    Examples
    --------
    >>> temp = np.array([11.28295376, 11.28295376])
    >>> em = np.array([4.78577516e+48, 4.78577516e+48])
    >>> rad_loss = calc_rad_loss(temp, em)
    >>> rad_loss
    array([  3.57994116e+26,   3.57994116e+26])
    """

    # Check inputs are correct
    temp = np.asanyarray(temp, dtype=np.float64)
    em = np.asanyarray(em, dtype=np.float64)

    # Initialize lists to hold model data of temperature - rad loss rate
    # relationship read in from csv file
    modeltemp = [] # modelled temperature is in log_10 sapce in units of MK
    model_loss_rate = []

    # Read data from csv file into lists, being sure to skip commented
    # lines begining with "#"
    with open(os.path.join(INSTR_FILES_PATH, "chianti_rad_loss.csv"),
              "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["temp_K"]))
            model_loss_rate.append(float(row["rad_loss_rate_per_em"]))
    modeltemp = np.asarray(modeltemp)
    model_loss_rate = np.asarray(model_loss_rate)
    # Ensure input values of flux ratio are within limits of model table
    if np.min(temp*1e6) < np.min(modeltemp) or
        np.max(temp*1e6) > np.max(modeltemp):
        raise ValueError("All values in temp must be within the range " +
                         "{0} - {1} MK.".format(np.min(modeltemp/1e6),
                                                np.max(modeltemp/1e6)))
    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modeltemp, model_loss_rate, s=0)
    rad_loss_rate = em * interpolate.splev(temp*1e6, spline, der=0)

    # If obstime keyword giving measurement times is set, calculate
    # radiative losses intergrated over time.
    if obstime is not None:
        dt = _time_intervals(obstime)
        # Check that times are in chronological order
        if np.min(dt) <= 0:
            raise ValueError("times in obstime must be in " +
                             "chronological order.")
        rad_loss_int = np.sum(rad_loss_rate*dt)
        rad_loss_out = {"temperature":temp, "em":em,
                        "rad_loss_rate":rad_loss_rate, "time": obstime,
                        "rad_loss_int":rad_loss_int}
    else:
        rad_loss_out = {"temperature":temp, "em":em,
                        "rad_loss_rate":rad_loss_rate}

    return rad_loss_out

def xray_luminosity(goeslc):
    """
    Calculates and adds solar X-ray luminosity to a GOESLightCurve.

    This function calculates the solar X-ray luminosity in the
    GOES wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based
    on the observed GOES fluxes.  The units of the results are erg/s.
    This is done by calling goes_lx().  This function assumes that the
    radiation is emitted isotropically, i.e. is distributed over a
    spherical surface area with a radius equal to the Sun-Earth
    distance.  Once the luminosity in each GOES passband is found,
    they are added to a copy of the original GOESLightCurve object as
    goeslc.data.luminosity_xrsa (for the 0.5-4 angstrom channel) and
    goeslc.data.luminosity_xrsb (for the 1-8 angstrom channel), where
    goeslc is the GOESLightCurve object.

    Parameters
    ----------
    goeslc : GOESLightCurve object

    Returns
    -------
    goeslc.data.luminosity_xrsa : pandas.core.series.Series
                                  Array of luminosity in the 0.5-4
                                  angstrom wavelength range [erg/s]
    goeslc.data.luminosity_xrsb : pandas.core.series.Series
                                  Array of luminosity in the 1-8
                                  angstrom wavelength range [erg/s]

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
                          xrsa   xrsb  luminosity_xrsa  luminosity_xrsb
    2014-01-01 00:00:00  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:02  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:04  7e-07  7e-06     1.903523e+24     1.903523e+25
    2014-01-01 00:00:06  7e-07  7e-06     1.903523e+24     1.903523e+25
    """
    # Check that input argument is of correct type
    if not isinstance(goeslc, sunpy.lightcurve.GOESLightCurve):
        raise TypeError("goeslc must be a GOESLightCurve object.")
    # Find temperature and emission measure with goes_chianti_tem
    lx_out = goes_lx(goeslc.data.xrsb, goeslc.data.xrsa,
                     date=str(goeslc.data.index[0]))
    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["luminosity_xrsa"] = lx_out["shortlum"]
    goeslc_new.data["luminosity_xrsb"] = lx_out["longlum"]

    return goeslc_new

def goes_lx(longflux, shortflux, obstime=None, date=None):
    """Calculates solar X-ray luminosity in GOES wavelength ranges.

    This function calculates the X-ray luminosity from the Sun in the
    GOES wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based
    on the observed GOES fluxes.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    longflux : ndarray or array-like which can be converted to float64
               type, such as an np.array, tuple, list.
               Array containing the observed GOES/XRS long channel flux
    shortflux : ndarray or array-like which can be converted to float64
                type, such as an np.array, tuple, list.
                Array containing the observed GOES/XRS short channel
                flux
    obstime : (optional) numpy ndarray, dtype=datetime64
              Measurement times corresponding to each long/short
              channel flux measurement.
    date : (optional) datetime object or valid date string
           Date at which measurements were taken.

    Returns
    -------
    longlum : numpy ndarray
              Array of luminosity in the long channel range
              (1-8 angstroms)
    shortlum : numpy ndarray
               Array of luminosity in the short channel range
               (0.5-4 angstroms)
    longlum_int : float
                  Long channel fluence, i.e. luminosity integrated
                  over time.
    shortlum_int : float
                   Short channel fluence, i.e. luminosity integrated
                   over time

    Notes
    -----
    This function calls goes_luminosity() to calculate luminosities.
    For more information on how this is done, see docstring of that
    function.

    Examples
    --------
    >>> longflux = np.array([7e-6,7e-6,7e-6,7e-6,7e-6,7e-6])
    >>> shortflux = np.array([7e-7,7e-7,7e-7,7e-7,7e-7,7e-7])
    >>> obstime = np.array(["2014-01-01 00:00:00",
                            "2014-01-01 00:00:02",
                            "2014-01-01 00:00:04"
                            "2014-01-01 00:00:06"
                            "2014-01-01 00:00:08"
                            "2014-01-01 00:00:10"]
                            dtype="datetime64[ms]")
    >>> lx_out = goes_lx(longflux, shortflux, obstime)
    >>> lx_out.longlum
    array([  1.98650769e+25,   1.98650769e+25,   1.98650769e+25,
             1.98650769e+25,   1.98650769e+25,   1.98650769e+25])
    >>> lx_out.shortlum
    array([  1.98650769e+24,   1.98650769e+24,   1.98650769e+24,
             1.98650769e+24,   1.98650769e+24,   1.98650769e+24])
    >>> lx_out.longlum_int
    2.0337865720138238e+26
    >>> lx_out.shortlum_int
    2.0337865720138235e+25

    """
    # Calculate X-ray luminosities
    longlum = _calc_xraylum(longflux, date=date)
    shortlum = _calc_xraylum(shortflux, date=date)

    # If obstime keyword giving measurement times is set, calculate
    # total energy radiated in the GOES bandpasses during the flare.
    if obstime is not None:
        dt = _time_intervals(obstime)
        # Check that times are in chronological order
        if np.min(dt) <= 0:
            raise ValueError("times in obstime must be in "
                             "chronological order.")
        longlum_int = np.sum(longlum*dt)
        shortlum_int = np.sum(shortlum*dt)
        lx_out = {"longflux":longflux, "shortflux":shortflux,
                  "time":obstime, "longlum":longlum, "shortlum":shortlum,
                  "longlum_int":longlum_int, "shortlum_int":shortlum_int,
                  "dt":dt}
    else:
        lx_out = {"longflux":longflux, "shortflux":shortflux,
                  "longlum":longlum, "shortlum":shortlum,}

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
    flux : ndarray or array-like which can be converted to float64
           type, such as an np.array, tuple, list.
           Array containing the observed solar flux
    date : (optional) datetime object or valid date string
           Used to calculate a more accurate Sun-Earth distance.

    Returns
    -------
    luminosity : numpy array
                 Array of luminosity

    Notes
    -----
    To convert from W/m**2 to erg/s:
    1 W = 1 J/s = 10**7 erg/s
    1 W/m**2 = 4*pi * AU**2 * 10**7 erg/s, where AU is the Sun-Earth
    distance in metres.

    Examples
    --------
    >>> flux = np.array([7e-6,7e-6])
    >>> luminosity = _calc_xraylum(flux, date="2014-04-21")
    >>> luminosity
    array([  1.98650769e+25,   1.98650769e+25])

    """
    # Ensure input is of correct type
    flux = np.asanyarray(flux, dtype=np.float64)
    if date is not None:
        date = parse_time(date) # Ensure date is of correct type
        return 4 * np.pi * (sun.constants.au.value * 
                            sun.sunearth_distance(t=date))**2 * 1e7 * flux
    else:
        return 4 * np.pi * (sun.constants.au.value)**2 * 1e7 * flux

def _check_download_file(filename, remotepath, localpath=os.path.curdir,
                         remotename=None, force_download=False):
    """
    Downloads a file from remotepath to localpath if it isn't there.

    This function checks whether a file with name filename exists in the
    location, localpath, on the user's local machine.  If it doesn't,
    it downloads the file from remotepath.

    Parameters
    ----------
    filename : string
               Name of file.
    remotepath : string
                 URL of the remote location from which filename can be
                 dowloaded.
    localpath : string
                Path of the directory in which filename should be
                stored.
                Default is current directory
    remotename : (optional) string
                 filename under which the file is stored remotely.
                 Default is same as filename.
    force_download : (optional) bool
                     If True, file will be downloaded whether or not
                     file already exists locally.

    Examples
    --------
    >>> pwd
    u'Users/user/Desktop/'
    >>> ls
    file.py
    >>> remotepath = "http://www.download_repository.com/downloads/"
    >>> _check_download_file("filename.txt", remotepath)
    >>> ls
    file.py    filename.txt

    """
    # Check if file already exists locally.  If not, try downloading it.
    if force_download or not os.path.isfile(os.path.join(localpath, filename)):
        # set local and remote file names be the same unless specified
        # by user.
        if type(remotename) is not str:
            remotename = filename
        try:
            # Check if the host server can be connected to.
            response = urllib2.urlopen(remotepath+remotename, timeout=1)
            # Try downloading file
            urllib.urlretrieve(os.path.join(remotepath, remotename),
                               os.path.join(localpath, filename))
            # Check if file has been downloaded.  If not, raise error.
            if not os.path.isfile(os.path.join(localpath, filename)):
                raise IOError(remotename + " was not downloaded from " +
                                remotepath + " .")
        except urllib2.URLError as e:
            # If the host server couldn't be connected to, raise Error.
            raise e(remotename + " could not be downloaded from\n" +
                    remotepath + " as connection could not be made.")
        except urllib2.HTTPError as e:
            # If file does not exist on server, print that info.
            if e.reason == "Not Found":
                raise e(remotename + " could be not be found at " + \
                  remotepath + ".")
            else:
                raise e

def _time_intervals(obstime):
    """
    Calculates time intervals between measurement times in seconds.

    This function calculates the time intervals between a series of
    measurement times for use in simple integration over time.
    Assume you have a series of times labelled t_1,...t_n.
    The start of the time bin for time t_i is defined as
    dt_i = (t_(i+1) - t_(i-1)) / 2
    i.e. from halfway between t_i and the previous time, t_(i-1), to
    halfway between t_i and the next time, t_(i+1).
    The time intervals for t_1 and t_n are special cases.  These are
    defined as
    dt_1 = (t_2 - t_1) / 2
    dt_(n-1) = (t_n - t_(n-1)) / 2
    In the case where only two time measurements are given, two
    "time intervals" are given which are equal to half the time
    difference between the two times given.

    Parameters
    ----------
    obstime : ndarray or array-like which can be converted to
              datetime64 type.
              Array containing the time measurements.

    Returns
    -------
    dt : numpy array, dtype=float
         Array of time intervals in [s]

    Examples
    --------
    >>> obstime = np.array(["2014-01-01 00:00:00",
                            "2014-01-01 00:00:02",
                            "2014-01-01 00:00:04"
                            "2014-01-01 00:00:06"
                            "2014-01-01 00:00:08"
                            "2014-01-01 00:00:10"],
                            dtype="datetime64[ms]")
    >>> dt = _time_intervals(obstime)
    >>> dt
    array([ 1.000,  2.000,  2.000,  2.000,  2.000,  1.000])

    """
    # check obstime is correct type and in units of milliseconds
    obstime = np.asarray(obstime, dtype="datetime64[ms]")
    # Ensure obstime has more than one element.  If so, calculate
    # difference between each time measurement.
    if len(obstime) < 2:
        raise IOError("obstime must have 2 or more elements.")
    elif len(obstime) == 2:
        dt = np.array((obstime[1]-obstime[0])/2)
        dt = np.append(dt, dt)
    else:
        dt = (obstime[2:]-obstime[:-2]) / 2
        dt = np.insert(dt, 0, (obstime[1]-obstime[0])/2)
        dt = np.append(dt, (obstime[-1]-obstime[-2])/2)
    # Finally, convert from [ms] to [s]
    dt = dt.astype(float) / 1e3
    return dt
