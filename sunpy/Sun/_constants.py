"""
Collection of solar physical constants. Most constants are in SI units.

The list is not meant to be comprehensive, but just a convenient list for everyday use.
"""

"""
Written: Steven Christe (7-May-2011)
Modified: 

physical constants: imported from Review of Particle Physics 2010 (page 102), NASA Sun Fact Sheet
(http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html), and 
Wikipedia (http://en.wikipedia.org/wiki/Sun)
Use at own risk:
TODO: References should be to published or standard sources (i.e. NOT websites)
"""

physical_constants = {}

# physical_constants[name] = (val, units, uncert)

physical_constants['mass'] = (1.9884e30, 'kg', -1)
physical_constants['radius'] = (6.9551e8, 'm', -1)
physical_constants['diameter'] = (6.9551e8*2.0, 'm', -1)
physical_constants['volume'] = (1.412e18, 'm^3', -1)
physical_constants['surface area'] = (6.0877e12, 'm^2', -1)
physical_constants['average density'] = ( 1.408e3, 'kg/m^3', -1)
physical_constants['center density'] = ( 1.622e5, 'kg/m^3', -1)
physical_constants['surface gravity'] = ( 274, 'kg/m^3', -1)
physical_constants['intensity'] = ( 2.009e7, '?', -1)
physical_constants['effective temperature'] = ( 5778, 'K', -1)
physical_constants['center temperature'] = ( 1.57e7, 'K', -1)
physical_constants['luminosity'] = ( 3.8427e26, 'J/s', -1)
physical_constants['absolute magnitude'] = ( 4.83, '', -1)
physical_constants['visual magnitude'] = ( -26.74, '', -1)
physical_constants['mass conversion rate'] = ( 4300e6, 'kg/s', -1)
physical_constants['mean energy production'] = ( 0.1937, 'J/kg', -1)
physical_constants['escape velocity'] = ( 617.6, 'km/s', -1)
physical_constants['ellipticity'] = ( 0.00005, '', -1)
physical_constants['GM'] = ( 132712e6, 'km^3/s^2', -1)
physical_constants['sunspot cycle'] = ( 11.4, 'years', -1)
physical_constants['metallicity'] = ( 0.0122, '', -1)