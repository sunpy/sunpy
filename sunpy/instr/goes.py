from __future__ import absolute_import
import sys

import numpy
import scipy.interpolate
import datetime
import dateutil
import csv

from sunpy.net import hek
from sunpy.time import parse_time

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
                          of time.  Must be of same length. [W/m^2].
    satellite : int
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    date : datetime object or str
           Date when observations made.  Important for correct
           calibration.  Default=today
    photospheric : bool
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

    Returns
    -------
    temp : numpy array
           Array of temperature values of same length as longflux and
           shortflux.  [MK]
    em : numpy array
         Array of volume emission measure values of same length as
         longflux and shortflux.  [10^49 cm^-3]

    Notes
    -----
    The temperature and volume emission measure are calculated here
    using the methods of White et al. (2005) who used the
    CHIANTI atomic physics database to model the response of the ratio
    of the short (0.5-4 angstrom) to long (1-8 angstrom) channels of the
    XRSs onboard various GOES satellites.  This method assumes an
    isothermal plasma, the ionisation equilibria of
    Mazzotta et al. (1998), and a constant density of 10^10 cm^-3.
    (See White et al. 2005 for justification of this last assumption.)
    This function is based on goes_chianti_tem.pro in SolarSoftWare
    written in IDL by Stephen White.

    Recent fluxes released to the public are scaled to be consistent
    with GOES-7.  In fact these recent fluxes are correct and so this
    correction must be removed before proceeding to use transfer
    functions.
    Email Rodney Viereck (NOAA) for more information.

    Measurements of short channel flux of less than 1e-10 W/m^2 or
    long channel flux less than 3e-8 W/m^2 are not considered good.
    Ratio values corresponding to suxh fluxes are set to 0.003.
         
    References
    ----------
    .. [1] White, S. M., Thomas, R. J., & Schwartz, R. A. 2005, Sol. Phys.,
       227, 231
    .. [2] Mazzotta, P., Mazzitelli, G., Colafrancesco, S., & Vittorio, N.
       1998, A&AS, 133, 339

    Examples
    --------
    >>> longflux = numpy.array([7e-6, 7e-6])
    >>> shortflux = numpy.array([7e-7, 7e-7])
    >>> temp, em = goes_chianti_tem(longflux, shortflux, satellite=15,
                                    date='2014-04-16', photospheric=False)
    >>> temp
    array([?????, ??????])
    >>> em
    array([?????, ??????])

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
    index = numpy.where(shortflux_corrected < 1e-10) or \
            numpy.where(longflux_corrected < 3e-8)
    fluxratio = shortflux_corrected / longflux_corrected
    fluxratio[index] = 0.003

    # FIND TEMPERATURE AND EMISSION MEASURE FROM FUNCTIONS BELOW
    temp = goes_get_chianti_temp(fluxratio, satellite=satellite,
                                 photospheric=photospheric)
    #em = goes_get_chianti_em(fluxratio, temp, satellite=satellite,
    #                         photospheric=photospheric)
    em = 1
    
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
    fluxratio : numpy ndarray
                Array containing the ratio of short channel to long
                channel GOES/XRS flux measurements.
    satellite : int
                Number of GOES satellite used to make observations.
                Important for correct calibration of data.
                Default=8
    photospheric : bool
                   States whether photospheric or coronal abundances
                   should be assumed.
                   Default=False, i.e. coronal abundances assumed.

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
    Mazzotta et al. (1998), and a constant density of 10^10 cm^-3.
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
    >>> fluxratio = numpy.array([10,10])
    >>> temp = goes_get_chianti_temp(fluxratio, satellite=15,
                                     photospheric=False)
    >>> temp
    array([?????, ???????])

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
    modeltemp = numpy.array(modeltemp)
    modelratio = numpy.array(modelratio)

    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = scipy.interpolate.splrep(modelratio, modeltemp, s=0)
    temp = 10.**scipy.interpolate.splev(fluxratio, spline, der=0)

    return temp
    
