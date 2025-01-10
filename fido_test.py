from sunpy.net import Fido, attrs as a
import astropy.units as u


result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Source.soho) 
print(result)

# Reproduces the bug where Fido.search crashes entirely even 
# if one of the client has an error using a FakeClient


