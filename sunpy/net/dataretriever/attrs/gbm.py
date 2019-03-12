from sunpy.net.attr import SimpleAttr

__all__ = ['Detector']

class Detector(SimpleAttr):
	"""
	Detector number for FERMI GBM
	"""


class Datatype(SimpleAttr):
	"""
	Data type of GBM - either CSPEC or CTIME
	"""