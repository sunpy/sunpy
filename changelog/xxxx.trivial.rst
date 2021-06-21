`sunpy.map.GenericMap.rsun_obs` no longer emits a warning if the metadata it
looks for is not present. Instead the standard photospheric radius is assumed
and a log message emitted at the 'info' level.
