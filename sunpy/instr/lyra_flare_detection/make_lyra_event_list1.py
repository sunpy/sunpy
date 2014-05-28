# Outline of code
#-------Prepare data-------
# 1.  Create LYRALightcurve object for given time interval from level 3
#     fits files (1-min averaged).
# 2.  Download LYTAF file from server.  Add keyword to allow user to
#     define local file.
# 3.  Remove data during UV occulations, offpoints, calibrations, and
#     large angle rotations.  (Also remove recovery times?)
# 4.  Download Ingolf's comb-like masks.  Add keyword so they can be
#     defined locally.
# 5.  Use the masks to remove comb-like features.
#------Find possible start times flares--------
# 6.  Take derivative of time series.
# 7.  Take where derivative is positive.
# 8.  Take where derivative is negative.
# 9.  Copy time array and shift by 4 elements.
# 10.  Subtract one from other and find where result = 4 minutes.
# 11. Define earlier time corresponding to previous step as preliminary
#     start times.
#-----Find end and peak times---------
# 12. Loop through start times.
# 13. Find next time where derivative in negative.
# 14. While end time criterion not reached, loop through following times
#     where derivative is negative and find when end criterion is
#     reached.
# 15. Define peak time as time of maximum between start and end time.
#-----Look for flares during decay phase------
# 16. If there are possible start times between peak and end time,
#     repeat steps 13-16.
def find_lyra_events(flux, time):
    """
    Finds events in a times series satisfying GOES event definitions.

    This function finds events in an input time series which satisfy
    the GOES event list definitions and returns the start, peak and end
    times.  The event list definitions are:
    Start time:
    1) There must be a continuous increase in 1-minute-averaged data
    over 4 minutes.
    2) The flux in the 4th minute must be at least 1.4 times the flux
    in the first minute.
    End time:
    1) The end time is when the flux falls to half-way between the peak
    and initial fluxes.

    Parameters
    ----------
    flux : numpy float64 array
           Contains flux/irradiance measurements
    time : sunpy datetime object or string array/list
           Contains measurement times corresponding to each element in
           flux.  Must be same length as flux.

    Returns
    -------
    event_list : 

    Refernces
    ---------
    http://www.swpc.noaa.gov/ftpdir/indices/events/README

    Examples
    --------                
    
    """
    # Define 
