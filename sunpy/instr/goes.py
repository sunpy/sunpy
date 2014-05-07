from __future__ import absolute_import
import sys

import numpy as np
import scipy.interpolate as interpolate
import datetime
import dateutil
import csv
import copy
from itertools import dropwhile

from sunpy.net import hek
from sunpy.time import parse_time
from sunpy.sun import sun
import sunpy.lightcurve
from sunpy.instr import exceptions

__all__ = ['get_goes_event_list']

def get_goes_event_list(trange,goes_class_filter=None):
    """A function to retrieve a list of flares detected by GOES within a given time range.

    Parameters
    ----------
    trange: a SunPy TimeRange object

    goes_class_filter: (optional) string
        a string specifying a minimum GOES class for inclusion in the list, e.g. 'M1', 'X2'.

    """
    
    
    #use HEK module to search for GOES events
    client=hek.HEKClient()
    event_type='FL'
    tstart=trange.start()
    tend=trange.end()

    #query the HEK for a list of events detected by the GOES instrument between tstart and tend (using a GOES-class filter)
    if goes_class_filter:
        result=client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),hek.attrs.FL.GOESCls > goes_class_filter,hek.attrs.OBS.Observatory == 'GOES')
    else:
        result=client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),hek.attrs.OBS.Observatory == 'GOES')

    #want to condense the results of the query into a more manageable dictionary
    #keep event data, start time, peak time, end time, GOES-class, location, active region source (as per GOES list standard)
    #make this into a list of dictionaries
    goes_event_list=[]

    for r in result:
        goes_event={}
        goes_event['event_date'] = parse_time(r['event_starttime']).date().strftime('%Y-%m-%d')
        goes_event['start_time'] =parse_time(r['event_starttime'])
        goes_event['peak_time'] = parse_time(r['event_peaktime'])
        goes_event['end_time'] = parse_time(r['event_endtime'])
        goes_event['goes_class'] = str(r['fl_goescls'])
        goes_event['goes_location'] = r['event_coord1'],r['event_coord2']
        goes_event['noaa_active_region'] = r['ar_noaanum']
        goes_event_list.append(goes_event)

    return goes_event_list

def temp_em(goeslc, photospheric=False):
    """
    Calculates and adds temperature and EM to a GOESLightCurve.

    Extended Summary
    ----------------
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
    photospheric : bool (optional)
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

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
    exceptions.check_goeslc(goeslc, varname="goeslc")

    # extract properties from GOESLightCurve object and change type to
    # that required by goes_chianti_em
    longflux = np.array(goeslc.data.xrsb)
    shortflux = np.array(goeslc.data.xrsa)
    satellite = int(goeslc.meta["TELESCOP"].split()[1])
    date = str(goeslc.data.index[0])

    # Find temperature and emission measure with goes_chianti_tem
    temp, em = goes_chianti_tem(longflux, shortflux, satellite=satellite,
                                date=date, photospheric=photospheric)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["temperature"] = temp
    goeslc_new.data["em"] = em

    return goeslc_new

def goes_chianti_tem(longflux, shortflux, satellite=8,
                     date=datetime.datetime.today(), photospheric=False):
    """
    Calculates temperature and emission measure from GOES/XRS data.

    Extended Summary
    ----------------
    This function calculates the isothermal temperature and volume
    emission measure of the solar soft X-ray emitting plasma observed by
    the GOES/XRS.  This is done using the observed flux ratio of the
    short (0.5-4 angstrom) to long (1-8 angstrom) channels.

    Parameters
    ----------
    longflux, shortflux : numpy ndarray
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
    photospheric : bool (optional)
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

    Returns
    -------
    temp : numpy ndarray
           Array of temperature values of same length as longflux and
           shortflux.  [MK]
    em : numpy ndarray
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
                                    date='2014-04-16', photospheric=False)
    >>> temp
    array([11.28295376, 11.28295376])
    >>> em
    array([  4.78577516e+48,   4.78577516e+48])

    """
    # CHECK INPUTS ARE CORRECT
    exceptions.check_float(longflux, varname="longflux") # Check longflux type
    exceptions.check_float(shortflux, varname="shortflux") # Check shortflux type
    satellite = exceptions.check_goessat(satellite) # Check satellite type
    date = exceptions.check_date(date) # Check date type
    # Check flux arrays are of same length.
    if len(longflux) != len(shortflux):
        raise ValueError("longflux and shortflux must have same number of " + \
                         "elements.")
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
    index = np.where(shortflux_corrected < 1e-10) or \
            np.where(longflux_corrected < 3e-8)
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio[index] = 0.003

    # FIND TEMPERATURE AND EMISSION MEASURE FROM FUNCTIONS BELOW
    temp = goes_get_chianti_temp(fluxratio, satellite=satellite,
                                 photospheric=photospheric)
    em = goes_get_chianti_em(longflux_corrected, temp, satellite=satellite,
                             photospheric=photospheric)
    return temp, em

def goes_get_chianti_temp(fluxratio, satellite=8, photospheric=False):
    """Calculates temperature from GOES flux ratio.

    Extended Summary
    ----------------
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
    fluxratio : numpy ndarray, dtype=float
                Array containing the ratio of short channel to long
                channel GOES/XRS flux measurements.
    satellite : int (optional)
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    photospheric : bool (optional)
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

    Returns
    -------
    temp : numpy ndarray
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
    >>> temp = goes_get_chianti_temp(fluxratio, satellite=15,
                                     photospheric=False)
    >>> temp
    array([11.28295376, 11.28295376])

    """

    # check inputs are correct
    exceptions.check_float(fluxratio, varname="fluxratio") # Check fluxratio type
    satellite = exceptions.check_goessat(satellite) # Check satellite type
    exceptions.check_photospheric(photospheric) # Check photospheric input

    # Initialize lists to hold model data of flux ratio - temperature
    # relationship read in from csv file
    modeltemp = [] # modelled temperature is in log_10 space in units of MK
    modelratio = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = "ratioGOES"+str(satellite)
    # Read data representing appropriate temperature--flux ratio
    # relationship depending on satellite number and assumed abundances.
    if photospheric is False:
        abund = "cor"
    else:
        abund = "pho"
    with open("goes_chianti_temp_"+abund+".csv", "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelratio.append(float(row[label]))
    modeltemp = np.array(modeltemp)
    modelratio = np.array(modelratio)

    # Ensure input values of flux ratio are within limits of model table
    if np.min(fluxratio) < np.min(modelratio) or \
      np.max(fluxratio) > np.max(modelratio):
        raise ValueError("For GOES " + str(satellite) + ", all values in " +
                         "fluxratio input must be within the range " +
                         str(np.min(modelratio)) + " - " +
                         str(np.max(modelratio)))

    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modelratio, modeltemp, s=0)
    temp = 10.**interpolate.splev(fluxratio, spline, der=0)

    return temp

def goes_get_chianti_em(longflux, temp, satellite=8, photospheric=False):
    """Calculates emission measure from GOES 1-8A flux and temperature.

    Extended Summary
    ----------------
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
    longflux : numpy ndarray, dtype=float
               Array containing the observed GOES/XRS long channel flux
    temp : numpy ndarray, dtype=float
           Array containing the GOES temperature
    satellite : int (optional)
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    photospheric : bool (optional)
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

    Returns
    -------
    em : numpy ndarray
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
    >>> em = goes_get_chianti_em(longflux, temp, satellite=15,
                                photospheric=False)
    >>> em
    array([  3.45200672e+48,   3.45200672e+48])

    """

    # Check inputs are correct
    exceptions.check_float(longflux, varname="longflux") # Check longflux input
    exceptions.check_float(temp, varname="temp") # Check temp input
    satellite = exceptions.check_goessat(satellite) # Check satellite type
    exceptions.check_photospheric(photospheric) # Check photospheric input
    # check input arrays are of same length
    if len(longflux) != len(temp):
        raise ValueError("longflux and temp must have same number of " + \
                         "elements.")

    # Initialize lists to hold model data of temperature - long channel
    # flux relationship read in from csv file.
    modeltemp = [] # modelled temperature is in log_10 sapce in units of MK
    modelflux = []
    # Determine name of column in csv file containing model ratio values
    # for relevant GOES satellite
    label = "longfluxGOES"+str(satellite)

    # Read data representing appropriate temperature--long flux
    # relationship depending on satellite number and assumed abundances.
    if photospheric is False:
        abund = "cor"
    else:
        abund = "pho"
    with open("goes_chianti_em_"+abund+".csv", "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelflux.append(float(row[label]))
    modeltemp = np.array(modeltemp)
    modelflux = np.array(modelflux)

    # Ensure input values of flux ratio are within limits of model table
    if np.min(np.log10(temp)) < np.min(modeltemp) or \
      np.max(np.log10(temp)) > np.max(modeltemp):
        raise ValueError("All values in temp must be within the range " +
                         str(np.min(10**modeltemp)) + " - " +
                         str(np.max(10**modeltemp)) + " MK.")

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(np.log10(temp), spline, der=0)
    em = longflux/denom * 1e55

    return em
