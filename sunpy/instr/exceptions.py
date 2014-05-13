import numpy as np
import datetime
import dateutil

import sunpy.lightcurve

def check_float(test, varname="This variable"):
    """
    Raises TypeError if test isn't a float type numpy array.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  (Printed if exception is raised.)
    """
    if type(varname) is not str:
        varname = "This variable"
    if type(test) is not np.ndarray or (test.dtype.type is not np.float64 and \
      test.dtype.type is not np.float32 and test.dtype.type is not np.float16):
        raise TypeError("{0} must be a numpy float array.".format(varname))

def check_goessat(test, varname="satellite"):
    """
    Raises Exception if test isn't a GOES satellite number.

    Raises a TypeError is test isn't an integer or can't be converted
    into a valid integer from a string.
    Raises a ValueError if test is an int less than 1.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  Default = 'satellite'
              (Printed if exception is raised.)

    Returns
    -------
    test : int
           Returned as original int if exceptions aren't raised, or a
           new int converted from input if input is a valid date string.

    """
    if type(varname) is not str:
        varname = "satellite"
    if type(test) is not int:
        if type(test) is str:
            try:
                test = int(test)
            except ValueError:
                raise TypeError("{0} must be an integer.".format(varname))
        else:
            raise TypeError("{0} must be an integer.".format(varname))
    if test < 1:
        raise ValueError("{0} must be the number (integer) of a "
                         "valid GOES satellite.".format(varname))
    return test

def check_photospheric(test, varname="photospheric"):
    """Raises TypeError if photospheric keyword isn't True or False.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  Default = 'photospheric'
              (Printed if exception is raised.)
    """
    if type(varname) is not str:
        varname = "photospheric"
    if type(test) is not bool:
        raise TypeError("{0} must be True or False.  \n"
                        "False: assume coronal abundances (default).  \n"
                        "True: assume photosperic abundances.".format(varname))

def check_date(test, varname="date"):
    """
    Raise TypeError if test isn't/can't be converted to datetime object.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  Default = 'date'
              (Printed if exception is raised.)

    Returns
    -------
    test : datetime object
           Returned as original datetime object if exceptions aren't
           raised, or a new datetime object converted from input if
           input is a valid date string.
           
    """
    if type(varname) is not str:
        varname = "date"
    if type(test) is not datetime.datetime:
        if type(test) is str:
            try: 
                test = dateutil.parser.parse(test)
            except TypeError:
                raise TypeError(
                    "{0} must be a datetime object.".format(varname))
        else:
            raise TypeError("{0} must be a datetime object.".format(varname))
    return test
    
def check_goeslc(test, varname="This variable"):
    """
    Raise TypeError if test is not a GOESLightCurve object.

    Parameters
    ----------
    test : variable to test
    varname : string, optional
              name of variable.  (Printed if exception is raised.)
              
    """
    if type(varname) is not str:
        varname = "This variable"
    if type(test) is not sunpy.lightcurve.sources.goes.GOESLightCurve:
        raise TypeError("{0} must be GOESLightCurve object.".format(varname))
