import numpy as np

import astropy.units as u

# TODO: Remove in astropy 5.1 (hopefully)
maxwell = u.def_unit(['Mx', 'Maxwell', 'maxwell'], 1e-8 * u.si.Wb, prefixes=True, doc="Maxwell: CGS unit for magnetic flux")
u.add_enabled_units([maxwell])

muSH = u.def_unit(
    "muSH",
    represents=((2 * np.pi * u.solRad**2)/1000000),
    prefixes=True,
    doc="Provide millionths of solar hemisphere units")
u.add_enabled_units([muSH])
