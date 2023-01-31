
from sunpy.net import Fido, attrs as a
#client = NOAAWeatherClient('differential-electrons-7-day.json')
results = Fido.search( a.Time("2016/1/1", "2016/1/2"), a.Instrument.goes)
print(results)
