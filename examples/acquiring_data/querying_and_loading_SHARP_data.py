"""
===============================
Querying and loading SHARP data
===============================

In this example we will demonstrate how to acquire [Spaceweather HMI Active Region Patch (SHARP)](http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm) data and load it into a `sunpy.map.Map`.
"""
import os

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###################################################################
# To search for SHARP data, we will need to query the [JSOC/Stanford](http://jsoc.stanford.edu/).
# We will use `Fido <sunpy.net.fido_factory.Fido>`
# and make use of the JSOC attributes that allow us to query the JSOC.

# You will need to replace this with your own registered email address:
# http://jsoc.stanford.edu/ajax/register_email.html
jsoc_email = os.environ["JSOC_EMAIL"]

result = Fido.search(a.Time("2023-02-01 21:00:00", "2023-02-01 22:30:00"),
                     a.Sample(1*u.hour),
                     a.jsoc.Series("hmi.sharp_cea_720s"),
                     a.jsoc.PrimeKey("HARPNUM", 7871),
                     a.jsoc.Notify(jsoc_email),
                     a.jsoc.Segment("Bp"))
print(result)

###################################################################
# Next, we can download the file

file = Fido.fetch(result)

###################################################################
# Now that we have the file, we can construct a `sunpy.map.Map` and plot it.

sharp_map = sunpy.map.Map(file)
sharp_map.plot()

plt.show()
