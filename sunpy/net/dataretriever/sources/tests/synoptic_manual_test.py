from sunpy.net import attrs as a, Fido
import astropy.units as u

req = a.Time("2020/01/01 00:00:00", "2020/01/05 00:00:00")
found = Fido.search(req, a.Instrument("AIASynoptic"), a.Sample(10 * u.minute))
print(found)
