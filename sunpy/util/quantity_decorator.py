# -*- coding: utf-8 -*-
import inspect

import astropy.units as u


class QuantityInput(object):
    """
    A decorator for a function that accepts some inputs as Quantity objects.
    
    This decorator attempts to convert the given Quantites to the units specified
    to the decorator, and fails nicely if a non-Quantity or a incompatible 
    unit was passed.
    
    Examples
    --------
    
    @QuantityInput(u.arcsec, u.arcsec, None)
    def myfunc(solarx, solary, someflag):
        pass
    """
    # __init__ is called when the function is parsed, it creates the decorator
    def __init__(self, **kwargs):
        self.equivalencies = kwargs.pop('equivalencies', [])
        self.f_kwargs = kwargs

    # __call__ is called when the wrapped function is called.
    def __call__(self, wrapped_function):
        
        #Define a new function to return in place of the wrapped one
        def wrapper(*func_args, **func_kwargs):
            
            args, varargs, keywords, defaults = inspect.getargspec(wrapped_function)

            for name, target_unit in self.f_kwargs.items():
                # We now parse the output of inspect to get the value of the
                # arg or the kwarg.
                if name in args:
                    loc = args.index(name)
                    # args includes names of kwargs which we want to skip here.
                    if loc < len(func_args):
                        arg = func_args[loc]

                elif name in func_kwargs or loc >= len(func_args):
                    arg = func_kwargs[name]
                
                else:
                    raise ValueError(
                        "Argument {0} is not a argument to wrapped function.".format(name))
                
                # Now we have the arg or the kwarg we check to see if it is 
                # convertable to the unit specified in the decorator.
                try:
                    arg.to(target_unit, equivalencies=self.equivalencies)

                # A UnitsError is raised if the conversion is not possible
                except u.UnitsError:
                    raise TypeError(
"Argument '{0}' to function '{1}' must be in units convertable to '{2}'.".format(
            name, wrapped_function.func_code.co_name, target_unit.to_string()))

                # AttributeError is raised if there is no `to` method.
                # i.e. not something that quacks like a Quantity.
                except AttributeError:
                    raise TypeError(
"Argument '{0}' to function '{1}' must be an astropy Quantity object".format(
                                     name, wrapped_function.func_code.co_name))

            return wrapped_function(*args, **func_kwargs)
        
        return wrapper
