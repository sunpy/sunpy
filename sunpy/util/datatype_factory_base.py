

class RegisteredFactoryBase(object):
    
    registry = dict()
    
    GenericWidgetType = None
    
    @classmethod
    def register(cls, WidgetType, dispatch_check_func=None, is_generic=False):
        
        if not is_generic:
            if dispatch_check_func is None:
                if hasattr(WidgetType, '_dispatch_check'):
                    cls.registry[WidgetType] = WidgetType._dispatch_check
            else:
                if callable(dispatch_check_func):
                    cls.registry[WidgetType] = dispatch_check_func
                else:
                    raise ValueError('dispatch_check_func must behave like a function')    
        else:
            cls.GenericWidgetType = WidgetType
    
    @classmethod
    def unregister(cls, WidgetType):
        cls.registry.pop(WidgetType)
                
    def __new__(cls, *args, **kwargs):
        
        if cls is RegisteredFactoryBase:
            
            WidgetType = None
            for key in cls.registry:
                
                if cls.registry[key](*args, **kwargs):
                    WidgetType = key
                    break
                    
            if WidgetType is not None:
                return WidgetType(*args, **kwargs)
            else:
                return cls.GenericWidgetType(*args, **kwargs)
        else:
            return super(RegisteredFactoryBase, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        raise NotImplementedError("RegistrationFactory __init__ should never be called.")