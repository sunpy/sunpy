import astropy.units as u

# TODO: Remove in astropy 5.1 (hopefully)
maxwell = u.def_unit(['Mx', 'Maxwell', 'maxwell'], 1e-8 * u.si.Wb, prefixes=True, doc="Maxwell: CGS unit for magnetic flux")
u.add_enabled_units([maxwell])
