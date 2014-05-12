def rad_loss_rate(goeslc):
    """
    Calculates and adds solar radiative loss rate to a GOESLightCurve.

    Extended Summary
    ----------------
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
    exceptions.check_goeslc(goeslc, varname="goeslc")

    # extract temperature and emission measure from GOESLightCurve
    # object and change type to that required by calc_rad_loss().
    # If GOESLightCurve object does not contain temperature and
    # emission measure, calculate using temp_em()
    try:
        temp = np.array(goeslc.data.temperature)
        em = np.array(goeslc.data.em)
    except AttributeError:
        goeslc_new = temp_em(goeslc)
        temp = np.array(goeslc_new.data.temperature)
        em = np.array(goeslc_new.data.em)
    else:
        goeslc_new = copy.deepcopy(goeslc)

    # Find radiative loss rate with calc_rad_loss()
    rad_loss_out = calc_rad_loss(temp, em)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new.data["rad_loss_rate"] = rad_loss_out["rad_loss_rate"]

    return goeslc_new

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

    # Check inputs are correct
    exceptions.check_float(temp, varname="temp") # Check temp type
    exceptions.check_float(em, varname="em") # Check em type

    # Initialize lists to hold model data of temperature - rad loss rate
    # relationship read in from csv file
    modeltemp = [] # modelled temperature is in log_10 sapce in units of MK
    model_loss_rate = []

    # Read data from csv file into lists, being sure to skip commented
    # lines begining with "#"
    with open(INSTR_FILES_PATH + "chianti_rad_loss.csv", "r") as csvfile:
        startline = dropwhile(lambda l: l.startswith("#"), csvfile)
        csvreader = csv.DictReader(startline, delimiter=";")
        for row in csvreader:
            modeltemp.append(float(row["temp_K"]))
            model_loss_rate.append(float(row["rad_loss_rate_per_em"]))
    modeltemp = np.array(modeltemp)
    model_loss_rate = np.array(model_loss_rate)
    # Ensure input values of flux ratio are within limits of model table
    if np.min(temp*1e6) < np.min(modeltemp) or \
      np.max(temp*1e6) > np.max(modeltemp):
        raise ValueError("All values in temp must be within the range " +
                         str(np.min(modeltemp/1e6)) + " - " +
                         str(np.max(modeltemp/1e6)) + " MK.")
    # Perform spline fit to model data to get temperatures for input
    # values of flux ratio
    spline = interpolate.splrep(modeltemp, model_loss_rate, s=0)
    rad_loss_rate = em * interpolate.splev(temp*1e6, spline, der=0)

    # If obstime keyword giving measurement times is set, calculate
    # radiative losses intergrated over time.
    if obstime is not None:
        dt = time_intervals(obstime)
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

    Extended Summary
    ----------------
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
    exceptions.check_goeslc(goeslc, varname="goeslc")

    # extract properties from GOESLightCurve object and change type to
    # that required by goes_chianti_em
    longflux = np.array(goeslc.data.xrsb)
    shortflux = np.array(goeslc.data.xrsa)
    date = str(goeslc.data.index[0])

    # Find temperature and emission measure with goes_chianti_tem
    lx_out = goes_lx(longflux, shortflux, date=date)

    # Enter results into new version of GOES LightCurve Object
    goeslc_new = copy.deepcopy(goeslc)
    goeslc_new.data["luminosity_xrsa"] = lx_out["shortlum"]
    goeslc_new.data["luminosity_xrsb"] = lx_out["longlum"]

    return goeslc_new

def goes_lx(longflux, shortflux, obstime=None, date=None):
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
    obstime : numpy ndarray, dtype=datetime64, optional
              Measurement times corresponding to each long/short
              channel flux measurement.
    date : datetime object or valid date strng, optional
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

    # Check inputs are of correct type
    exceptions.check_float(longflux, varname="longflux") # Check longflux type
    exceptions.check_float(shortflux, varname="shortflux") # Check shortflux type

    # Calculate X-ray luminosities
    longlum = calc_xraylum(longflux, date=date)
    shortlum = calc_xraylum(shortflux, date=date)

    # If obstime keyword giving measurement times is set, calculate
    # total energy radiated in the GOES bandpasses during the flare.
    if obstime is not None:
        dt = time_intervals(obstime)
        # Check that times are in chronological order
        if np.min(dt) <= 0:
            raise ValueError("times in obstime must be in " + \
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

def calc_xraylum(flux, date=None):
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
    date : datetime object or valid date string, optional
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
    >>> luminosity = goes_luminosity(flux, date="2014-04-21")
    >>> luminosity
    array([  1.98650769e+25,   1.98650769e+25])

    """
    exceptions.check_float(flux)
    if date is not None:
        date = exceptions.check_date(date)
        return 4 * np.pi * (sun.constants.au.value * 
                            sun.sunearth_distance(t=date))**2 * 1e7 * flux
    else:
        return 4 * np.pi * (sun.constants.au.value)**2 * 1e7 * flux

def time_intervals(obstime):
    """
    Calculates time intervals between measurement times in seconds.

    Extended Summary
    ----------------
    This function calculates the time intervals between a series of
    measurement times for use in siple integration over time.
    Assume you have a series of times labelled t_1,...t_n.
    The start of the time bin for time t_i is defined as
    dt_i = (t_(i+1) - t_(i-1)) / 2
    i.e. from halfway between t_i and the previous time, t_(i-1), to
    halfway between t_i and the next time, t_(i+1).
    The time intervals for t_1 and t_n are special cases.  These are
    defined as
    dt_1 = (t_2 - t_1) / 2
    dt_(n-1) = (t_n - t_(n-1)) / 2

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
    # check obstime is correct type and greater than min required length
    exceptions.check_datetime64(obstime, varname="obstime")
    if len(obstime) < 3:
        raise ValueError("obstime must have 3 or more elements")
    obstime = obstime.astype("datetime64[ms]")  # convert to units of ms
    dt = (obstime[2:]-obstime[:-2]) / 2
    dt = np.insert(dt, 0, (obstime[1]-obstime[0])/2)
    dt = np.append(dt, (obstime[-1]-obstime[-2])/2)
    dt = dt.astype(float) / 1e3 # convert from [ms] to [s]
    return dt

def check_datetime64(test, varname="This variable"):
    """Raise Exception if test isn't numpy array of dtype datetime64.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  (Printed if exception is raised.)
    """
    if type(varname) is not str:
        varname = "This variable"
    if type(test) is not np.ndarray or test.dtype.type is not np.datetime64:
        raise TypeError(varname + " must be a numpy array of type datetime64.")
