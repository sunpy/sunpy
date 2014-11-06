# -*- coding: utf-8 -*-

import astropy.units as u

def quantity_input(*f_args, **f_kwargs):
    """
    A decorator for a function that accepts some inputs as Quantity objects.
    
    This decorator attempts to convert the given Quantites to the units specified
    to the decorator, and fails nicely if a non-Quantity or a incompatible 
    unit was passed.
    
    Examples
    --------
    
    @quantity_input(u.arcsec, u.arcsec, None)
    def myfunc(solarx, solary, someflag):
        pass
    """
    equivalencies = f_kwargs.pop('equivalencies', [])
    def check_quantities(f):
        # Number of args in decorator must equal number of args in function
        if f.func_defaults:
            num_kwargs = len(f.func_defaults)
        else:
            num_kwargs = 0
        num_args = f.func_code.co_argcount - num_kwargs
        if len(f_args) != num_args:
            raise ValueError("Number of decorator arguments does not equal number of function arguments")

        def new_f(*args, **kwds):
            # Check args, number of args in decorator must equal number of args in function
            args = list(args)
            for i, (arg, f_arg) in enumerate(zip(args, f_args)):
                if f_arg is not None:
                    try:
                        arg.to(f_arg, equivalencies=equivalencies)
                    except u.UnitsError:
                        raise TypeError("Argument '{}' to function '{}' must be in units convertable to '{}'.".format(
                                        f.func_code.co_varnames[i], f.func_code.co_name, f_arg.to_string()))
                    except AttributeError:
                        raise TypeError("Argument '{}' to function '{}' must be an astropy Quantity object".format(
                                        f.func_code.co_varnames[i], f.func_code.co_name))

            # Check kwargs, only kwargs specified in the decorator are modified
            for j, (kwarg, value) in enumerate(f_kwargs.items()):
                if kwarg in kwds:
                    try:
                        kwds[kwarg].to(value, equivalencies=equivalencies)
                    except u.UnitsError:
                        raise TypeError("Keyword argument '{}' to function '{}' must be in units convertable to '{}'.".format(
                                        kwarg, f.func_code.co_name, value.to_string()))
                    except AttributeError:
                        raise TypeError("Argument '{}' to function '{}' must be an astropy Quantity object".format(
                                        f.func_code.co_varnames[j+num_args], f.func_code.co_name))

            return f(*args, **kwds)


        new_f.func_name = f.func_name
        return new_f

    return check_quantities
