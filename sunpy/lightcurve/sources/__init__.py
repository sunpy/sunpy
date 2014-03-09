"""
Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.LightCurve` class.
"""

from .. lightcurve_factory import LightCurve

from . lyra import LYRALightCurve
LightCurve.register(LYRALightCurve, LYRALightCurve._is_datasource_for)

from . goes import GOESLightCurve
LightCurve.register(GOESLightCurve, GOESLightCurve._is_datasource_for)

from . eve import EVELightCurve
LightCurve.register(EVELightCurve, EVELightCurve._is_datasource_for)

from . norh import NoRHLightCurve
LightCurve.register(NoRHLightCurve, NoRHLightCurve._is_datasource_for)

from . logical import LogicalLightCurve
LightCurve.register(LogicalLightCurve, LogicalLightCurve._is_datasource_for)