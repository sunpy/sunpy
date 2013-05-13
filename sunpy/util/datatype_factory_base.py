

class RegisteredFactoryBase(object):
	
	registry = dict()
	
	GenericWidgetType = None
	
	@classmethod
	def register(WidgetType, dispatch_check_func=None, is_generic=False):
		
		if not is_generic:
			if type_match_func is None:
				if hasattr(WidgetType, '_dispatch_check'):
					self.registry[WidgetType] = WidgetType._dispatch_check
			else:
				if callable(dispatch_check_func):
					self.registry[WidgetType] = dispatch_check_func
				else:
					raise ValueError('dispatch_check_func must behave like a function')	
		else:
			self.GenericWidgetType = WidgetType
	
	@classmethod
	def unregister(WidgetType):
		self.registry.pop(WidgetType)
				
	def __new__(cls, *args, **kwargs):
		
		if cls is RegistrationFactoryBase:
			
			WidgetType = None
			for key in self.registry:
				
				if self.registry[key](*args, **kwargs):
					WidgetType = key
					break
					
			if WidgetType is not None:
				return WidgetType(*args, **kwargs)
			else:
				return self.GenericWidgetType(*args, **kwargs)
		else:
			return super(RegistrationFactoryBase, cls).__new__(cls, *args, **kwargs)

	def __init__(self, *args, **kwargs):
		raise NotImplementedError("RegistrationFactory __init__ should never be called.")