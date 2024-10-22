from sunpy.net import Fido
from sunpy.net import attrs as a
import astropy.units as u

# Example usage
time_range = a.Time("2013-10-20 00:00:00", "2023-10-21 00:00:00")
instrument = a.Instrument("AIA")
wavelength = a.Wavelength(193 * u.angstrom)
sample = a.Sample(30*u.day)  #TBD
results = Fido.search(time_range, instrument, wavelength, a.Level("synoptic"), sample)   #, a.Level("synoptic"))  #, sample)
print(results)
