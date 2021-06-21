`sunpy.map.GenericMap.rsun_meters` now uses `sunpy.map.GenericMap.rsun_obs`
as a fallback to calculate the assumed radius of emission if RSUN_REF metadata
isn't present but metadata for `~sunpy.map.GenericMap.rsun_obs` is.
