class RegisteredFactoryBase(object):
    """
    Base class for factories with registerable widgets.  This allows classes 
    developed _outside_ of SunPy to be accessed via the same factories as 
    SunPy's built-in data types.
    
    Attributes
    ----------
    registry: dict
        The dictionary containing the mapping of WidgetType's to functions that
        match the arguments to that WidgetType.
    DefaultWidgetType: type
        Fall-back type for the event that no WidgetTypes match the specified 
        arguments.
    """
    
    registry = dict()
    
    DefaultWidgetType = None
    
    @classmethod
    def register(cls, WidgetType, dispatch_check_func=None, is_default=False):
        """ 
        Register a WidgetType with the factory.
        
        Parameters
        ----------
        WidgetType: type
            The type to register.
        dispatch_check_func: function-like, optional
            Function-like which returns True if arguments match WidgetType, 
            otherwise False.
        is_default: boolean, optional
            Allows specification of cls.DefaultWidgetType
            
        """ 
           
        if not is_default:
            # If no separate function is specified to match, check for a default
            # attribute.
            if dispatch_check_func is None:
                if hasattr(WidgetType, '_dispatch_check'):
                    if not callable(WidgetType._dispatch_check):
                        raise ValueError("""WidgetType._dispatch_check must
                                            behave like a function""")    
                    cls.registry[WidgetType] = WidgetType._dispatch_check
                else:
                    raise AttributeError("""dispatch_check_func must be specified 
                                          or WidgetType must have 
                                          _dispatch_check attribute.""")
            else:
                if callable(dispatch_check_func):
                    cls.registry[WidgetType] = dispatch_check_func
                else:
                    raise ValueError("""dispatch_check_func must behave like a 
                                      function""")    
        else:
            cls.DefaultWidgetType = WidgetType
    
    @classmethod
    def unregister(cls, WidgetType):
        """ Remove a type from the factory. """
        cls.registry.pop(WidgetType)
                
    def __new__(cls, *args, **kwargs):
        
        """
        Factory-construction method.  Takes arbitrary arguments and keyword 
        arguments and passes them to a sequence of pre-registered types to
        determine which is the correct Widget to build.
        """
        
        if cls is RegisteredFactoryBase:
            
            WidgetType = None
            
            # Loop over each registered type and check to see if WidgetType
            # matches the arguments.  If it does, use that type.
            for key in cls.registry:
                if cls.registry[key](*args, **kwargs):
                    WidgetType = key
                    break
            
            # If a type matches, instantiate it, otherwise fall back to the
            # default type.
            if WidgetType is not None:
                return WidgetType(*args, **kwargs)
            else:
                return cls.DefaultWidgetType(*args, **kwargs)
        else:
            return super(RegisteredFactoryBase, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        raise NotImplementedError("{0} __init__ should never be called.".format(self.__class__.__name__))