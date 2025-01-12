from sunpy.net import Fido
from sunpy.net import attrs as a

result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lasco)
print(result)

print("Errors dict", result.errors)
