from __future__ import absolute_import
from __future__ import division

import os.path
import datetime
import csv
import copy
import socket
from itertools import dropwhile

import numpy as np
import scipy.interpolate as interpolate

from sunpy.net import hek
from sunpy.time import parse_time
from sunpy import config
from sunpy import lightcurve
from sunpy.util.net import check_download_file

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

def temp_em(goeslc, abundances="coronal", download=False, download_dir=DATA_PATH):
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

    abundances : (optional) string equalling either 'coronal' or 'photospheric'.
        States whether photospheric or coronal abundances should be assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES temperature and emission measure data files are downloaded.
        It is important to do this if a new version of the files has been
        generated due to a new CHIANTI version being released or the launch of
        new GOES satellites since these files were originally downloaded.
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
    longflux, shortflux : ndarray or array-like which can be converted to float64 type, such as an np.array, tuple, list.
        Arrays containing the long and short GOES/XRS flux measurements
        respectively as a function of time.  Must be of same length. [W/m**2].

    satellite : int (optional)
        Number of GOES satellite used to make observations, important for
        correct calibration of data.
        Default=8

    date : datetime object or str
         Date when observations made.  Important for correctcalibration.
         Default=today

    abundances : (optional) string equalling either 'coronal' or 'photospheric'.
        States whether photospheric or coronal abundances should be assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES temperature and emission measure data files are
        downloaded.  It is important to do this if a new version of the files
        has been generated due to a new CHIANTI version being released or the
        launch of new GOES satellites since these files were originally downloaded.
        Default=False

    download_dir : (optional) string
        The directory to download the GOES temperature and emission measure
        data files to.
        Default=SunPy default download directory

    Returns
    -------
    temp : numpy array
        Array of temperature values of same length as longflux and shortflux. [MK]

    em : numpy array
        Array of volume emission measure values of same length as longflux
        and shortflux.  [10**49 cm**-3]

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
    fluxratio : ndarray or array-like which can be converted to float64 type, such as an np.array, tuple, list.
        Array containing the ratio of short channel to long channel GOES/XRS
        flux measurements.

    satellite : int (optional)
        Number of GOES satellite used to make observations. Important for
        correct calibration of data.
        Default=8

    abundances : (optional) string equalling either 'coronal' or 'photospheric'.
        States whether photospheric or coronal abundances should be assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES temperature data files are downloaded.
        It is important to do this if a new version of the files has been
        generated due to a new CHIANTI version being released or the launch
        of new GOES satellites since these files were originally downloaded.
        Default=False

    download_dir : (optional) string
        The directory to download the GOES temperature and emission measure
        data files to.
        Default=SunPy default download directory

    Returns
    -------
    temp : numpy array
        Array of temperature values of same length as longflux and shortflux. [MK]

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
    check_download_file(FILE_TEMP_COR, GOES_REMOTE_PATH, download_dir,
                        replace=download)
    check_download_file(FILE_TEMP_PHO, GOES_REMOTE_PATH, download_dir,
                        replace=download)

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
    longflux : ndarray or array-like which can be converted to float64 type, such as an np.array, tuple, list.
        Array containing the observed GOES/XRS long channel flux

    temp : ndarray or array-like which can be converted to float64 type, such as an np.array, tuple, list.
        Array containing the GOES temperature

    satellite : int (optional)
        Number of GOES satellite used to make observations.
        Important for correct calibration of data.
        Default=8

    abundances : (optional) string equalling either 'coronal' or 'photospheric'.
        States whether photospheric or coronal abundances should be assumed.
        Default='coronal'

    download : (optional) bool
        If True, the GOES emission measure data files are downloaded.
        It is important to do this if a new version of the files has been
        generated due to a new CHIANTI version being released or the launch of
        new GOES satellites since these files were originally downloaded.
        Default=False

    download_dir : (optional) string
        The directory to download the GOES temperature and emission measure
        data files to.
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
    check_download_file(FILE_EM_COR, GOES_REMOTE_PATH, download_dir,
                        replace=download)
    check_download_file(FILE_EM_PHO, GOES_REMOTE_PATH, download_dir,
                        replace=download)

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
      np.max(np.log10(temp)) > np.max(modeltemp) or \
      np.isnan(np.min(np.log10(temp))):
        raise ValueError("All values in temp must be within the range "
                         "{0} - {1} MK.".format(np.min(10**modeltemp),
                                                np.max(10**modeltemp)))

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(np.log10(temp), spline, der=0)
    em = longflux/denom * 1e55

    return em
