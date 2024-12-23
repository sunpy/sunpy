Fixes a bug where `~sunpy.map.sources.HMIMap` returned a wavelength without a unit because ``WAVEUNIT``
is not in the header and cannot be parsed from any other part of the metadata. If it cannot be found,
it now defaults to Angstrom.
