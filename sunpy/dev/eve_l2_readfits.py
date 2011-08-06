
# http://lasp.colorado.edu/eve/data_access/evewebdataproducts/level2/EVE_L2_V2_README.pdf

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

import datetime
import pyfits

def eve_l2_readfits(filename):
    """Read and interpret and EVE L2 fits file (each Level 2 data file spans one hour).
    Returns a dictionary.
    Documentation describing the structure of these files can be found
    
    http://lasp.colorado.edu/eve/data_access/evewebdataproducts/level2/EVE_L2_V2_README.pdf
    """
    
    filename = '/Users/schriste/Dropbox/EVL_L2_2011030_00_002_01.fit.fit.fit.gz'
    fits = pyfits.open(filename)
    
    science_data = fits[5].data
    
    # TODO: Should convert tai to datetimes with anytim once implemented and check with times
    #tai = science_data.field(0)
    yyyydoy = science_data.field(1)
    seconds_of_day = science_data.field(2)
    
    # the following checked with IDL anytim
    temp = zip(seconds_of_day, yyyydoy)
    years = [datetime.datetime(int(str(s[1])[0:4]), 1, 1) for s in temp]
    dt = [datetime.timedelta(int(str(s[1])[4:7]) - 1, s[0], 0) for s in temp]
    times = [s[0] + s[1] for s in zip(years, dt)]
    
    flags = science_data.field(3)   # 0 = good
    sc_flags = science_data.field(4)
    line_irradiance = science_data.field(5)
    line_precision = science_data.field(6)
    line_accuracy = science_data.field(7)
    band_irradiance = science_data.field(8)
    band_precision = science_data.field(9)
    band_accuracy = science_data.field(10)
    diode_irradiance = science_data.field(11)
    diode_stdev = science_data.field(12)
    diode_precision = science_data.field(13)
    quad_fraction = science_data.field(14)
    quad_stdev = science_data.field(15)
    quad_precision = science_data.field(16)
    
    # gives the units for everything above
    science_data_units = fits[6].data
    
    linesmeta = fits[1].data
    
    wave_center = linesmeta.field(0)
    wave_min = linesmeta.field(1)
    wave_max = linesmeta.field(2)
    logt = linesmeta.field(3)
    line_name = linesmeta.field(4)
    line_type = linesmeta.field(5)
    line_blends = linesmeta.field(6)
    
    bandsmeta = fits[2].data
    
    band_name = bandsmeta.field(0)
    band_type = bandsmeta.field(1)
    
    diodemeta = fits[3].data
    
    diode_name = diodemeta.field(0)
    diode_type = diodemeta.field(1)
    
    quadmeta = fits[4].data
    
    quad_name = quadmeta.field(0)
    quad_type = quadmeta.field(1)
    
    lines = []
    for i in range(len(line_name) - 1):
        line_data = {'name': line_name[i], 'type': line_type[i], 'blends': line_blends[i], 'wave_center': wave_center[i], 'wave_min': wave_min[i], 'wave_max': wave_max[i], 'log_temp': logt[i], 'irradiance': line_irradiance[i], 'precision': line_precision[i], 'accuracy': line_accuracy[i]}
        lines.append(line_data)
    
    bands = []
    for i in range(len(band_name) - 1):
        band_data = {'name': band_name[i], 'type': band_type[i], 'irradiance': band_irradiance[i], 'precision': band_precision[i], 'accuracy': band_accuracy[i]}
        bands.append(band_data)
        
    diodes = []
    for i in range(len(diode_name) - 1):
        diode_data = {'name': diode_name[i], 'type': diode_type[i], 'irradiance': diode_irradiance[i], 'stdev': diode_stdev[i], 'precision': diode_precision[i]}
        diodes.append(diode_data)
    
    quads = []
    for i in range(len(quad_name) - 1):
        diode_data = {'name': quad_name[i], 'type': quad_type[i], 'fraction': quad_fraction[i], 'stdev': quad_stdev[i], 'precision': quad_precision[i]}
        quads.append(diode_data)
        
    eve_l2 = {'times': times, 'flags': flags, 'sc_flags': sc_flags, 'units': science_data_units, 'lines': lines, 'bands': bands, 'diodes': diodes, 'quads': quads}
    
    return eve_l2

