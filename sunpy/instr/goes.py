from __future__ import absolute_import
import sys

import numpy as np
import scipy.interpolate as interpolate
import datetime
import dateutil
import csv
from itertools import dropwhile

from sunpy.net import hek
from sunpy.time import parse_time
from sunpy.sun import sun

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

    # CHECK INPUTS ARE OF CORRECT TYPE
    # Check flux arrays are of same length.
    if len(longflux) != len(shortflux):
        sys.exit("longflux and shortflux are not of same length.  Try again " +
                 "with longflux an shortflux of same length.")
    # Check satellite is an int and greater than zero
    while type(satellite) is not int:
        try:
            satellite = int(satellite)
        except ValueError:
            print "Attention: satellite must be an integer.  Try again."
            satellite = int(raw_input("Enter the number of the GOES " +
                                      "satellite as an integer: "))
        if type(satellite) is int and satellite < 1:
            print "Attention: satellite must be greater than 0.  Try again."
            satellite = raw_input("Enter GOES satellite number as integer "
                                  "greater than 0: ")
    # Ensure date is a datetime object or a string.
    # If date is a string, convert it to a datetime object
    while type(date) is not datetime.datetime and type(date) is not str:
        print "Attention: date must be a datetime object or string.  Try again."
        date = raw_input("Enter date as datetime object or string: ")
    if type(date) is str:
        date = dateutil.parser.parse(date)
    # Ensure photospheric is a Boolean value
    while type(photospheric) is not bool:
        print "Attention: photospheric must be a boolean value."
        photospheric = raw_input("Do you want to assume photospheric or " +
                                 "coronal abundances? Enter True for " +
                                 "photospheric and False for coronal.")
    
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
    # Initialize lists to hold model data of flux ratio - temperature
    # relationship read in from csv file
    modeltemp = [ ] # modelled temperature is in log_10 sapce in units of MK
    modelratio = [ ]
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
        csvreader = csv.DictReader(csvfile, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelratio.append(float(row[label]))
    modeltemp = np.array(modeltemp)
    modelratio = np.array(modelratio)

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

    # Initialize lists to hold model data of temperature - long channel
    # flux relationship read in from csv file.
    modeltemp = [ ] # modelled temperature is in log_10 sapce in units of MK
    modelflux = [ ]
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
        csvreader = csv.DictReader(csvfile, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["log10temp_MK"]))
            modelflux.append(float(row[label]))
    modeltemp = np.array(modeltemp)
    modelflux = np.array(modelflux)

    # Perform spline fit to model data
    spline = interpolate.splrep(modeltemp, modelflux, s=0)
    denom = interpolate.splev(np.log10(temp), spline, der=0)
    em = longflux/denom * 1e55

    return em

def calc_rad_loss(temp, em, obstime=None):
    """
    Finds radiative loss rate of solar SXR plasma over all wavelengths.

    Extended Summary
    ----------------
    This function calculates the luminosity of solar coronal soft
    X-ray-emitting plasma across all wavelengths given an isothermal
    temperature and emission measure.  The units of the results are
    erg/s.  This function is based on calc_rad_loss.pro in SSW IDL.
    In addition, if obstime keyword is set, giving the times to which
    the temperature and emission measure values correspond, the
    radiated losses integrated over time are also calculated.

    Parameters
    ----------
    temp : numpy ndarray, dtype=float, units=[MK]
           Array containing the temperature of the coronal plasma at
           different times.
    em : numpy ndarray, dtype=float, units=[cm**-3]
         Array containing the emission measure of the coronal plasma
         at the same times corresponding to the temperatures in temp.
         Must be same length as temp
    obstime : numpy array, dtype=datetime64, optional
              array of measurement times to which temperature and
              emission measure values correspond.  Must be same length
              as temp and em.  If this keyword is set, the integrated
              radiated energy is calculated.
    
    Returns
    -------
    radlossrate : numpy ndarray, dtype=float
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

    # Initialize lists to hold model data of temperature - rad loss rate
    # relationship read in from csv file
    model_temp = [ ] # modelled temperature is in log_10 sapce in units of MK
    model_loss_rate = [ ]

    # Read data from csv file into lists, being sure to skip commented
    # lines begining with "#"
    with open("chianti_rad_loss.csv", "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            model_temp.append(float(row["temp_K"]))
            model_loss_rate.append(float(row["rad_loss_rate_per_em"]))
    model_temp = np.array(model_temp)
    model_loss_rate = np.array(model_loss_rate)
    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(model_temp, model_loss_rate, s=0)
    rad_loss_rate = em * interpolate.splev(temp*1e6, spline, der=0)

    # If obstime keyword giving measurement times is set, calculate
    # radiative losses intergrated over time.
    if obstime is not None:
        try:
            timetype = obstime.dtype.type
        except AttributeError:
            timetype = None
        if timetype is np.datetime64:
            dt = time_intervals(obstime)
            rad_loss_int = np.sum(rad_loss_rate*dt)
            rad_loss_out = {"rad_loss_rate":rad_loss_rate,
                            "time": obstime,
                            "rad_loss_int":rad_loss_int}
        else:
             print "Warning: obstime must be of type datetime64.  \n" + \
               "         Integrated radiated losses will not be calculated."
             rad_loss_out = {"rad_loss_rate":rad_loss_rate}
    else:
        rad_loss_out = {"rad_loss_rate":rad_loss_rate}

    return rad_loss_out

def goes_lx(longflux, shortflux, obstime, date=None):
    """Calculates solar X-ray luminosity in GOES wavelength ranges.

    Extended Summary
    ----------------
    This function calculates the X-ray luminosity from the Sun in the
    GOES wavelength ranges (1-8 angstroms and 0.5-4 angstroms) based
    on the observed GOES fluxes.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    longflux : numpy ndarray, dtype=float
               Array containing the observed GOES/XRS long channel flux
    shortflux : numpy ndarray, dtype=float
                Array containing the observed GOES/XRS short channel
                flux
    obstime : numpy ndarray, dtype=datetime64
              Measurement times corresponding to each long/short
              channel flux measurement.

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
    longlum = goes_luminosity(longflux, date=date)
    shortlum = goes_luminosity(shortflux, date=date)

    # Next calculate the total energy radiated in the GOES bandpasses
    # during the flare.
    # First determine time intervals over which to intergrate each flux
    # measurement using time_intervals() function (below).
    dt = time_intervals(obstime)
    # Calculate integrated X-ray radiative losses over time duration.
    longlum_int = np.sum(longlum*dt)
    shortlum_int = np.sum(shortlum*dt)
    # put data together in a dictionary for output
    lx_out = {"longflux":longflux, "shortflux":shortflux, "time":obstime,
              "longlum":longlum, "shortlum":shortlum,
              "longlum_int":longlum_int, "shortlum_int":shortlum_int, "dt":dt}
    
    return lx_out

def goes_luminosity(flux, date=None):
    """
    Calculates solar luminosity based on observed flux observed at 1AU.

    Extended Summary
    ----------------
    This function calculates the luminosity from the Sun based
    on observed flux in W/m**2.  The units of the results are erg/s.
    The calculation is made by simply assuming that the radiation is
    emitted isotropically, i.e. is distributed over a spherical
    surface area with a radius equal to the Sun-Earth distance.

    Parameters
    ----------
    flux : numpy ndarray
           Array containing the observed solar flux

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
    >>> luminosity = goes_luminosity(flux, date="2014-04-21")
    >>> luminosity
    array([  1.98650769e+25,   1.98650769e+25])

    """
    return 4 * np.pi * \
      (sun.constants.au.value * sun.sunearth_distance(t=date))**2 * 1e7 * flux

def time_intervals(obstime):
    """
    Calculates time intervals between measurement times in seconds.

    Extended Summary
    ----------------
    This function calculates the time intervals between a series of
    measurement times for use in siple integration over time.
    Assume you have a series of times labelled t_1,...t_n.
    The start of the time bin for time t_i is defined as
    dt_i = ( t_(i+1) - t_(i-1) ) / 2
    i.e. from halfway between t_i and the previous time, t_(i-1), to
    halfway between t_i and the next time, t_(i+1).
    The time intervals for t_1 and t_n are special cases.  These are
    defined as
    dt_1 = ( t_2 - t_1 ) / 2
    dt_(n-1) = ( t_n - t_(n-1) ) / 2
    

    Parameters
    ----------
    obstime : numpy ndarray, dtype=datetime64
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
    >>> dt = time_intervals(obstime)
    >>> dt
    array([ 1.000,  2.000,  2.000,  2.000,  2.000,  1.000])

    """

    obstime = obstime.astype("datetime64[ms]")  # convert to units of ms
    dt = (obstime[2:]-obstime[:-2]) / 2
    dt = np.insert(dt, 0, (obstime[1]-obstime[0])/2)
    dt = np.append(dt, (obstime[-1]-obstime[-2])/2)
    dt = dt.astype(float) / 1e3 # convert from [ms] to [s]
    return dt
