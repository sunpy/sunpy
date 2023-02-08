from solarmach import SolarMACH
import datetime

# necessary options
body_list = ['STEREO-A', 'Earth', 'Mars']
vsw_list = [400, 400, 400]  
time = datetime.datetime.now()
date = str(time)

# optional parameters
coord_sys = 'Stonyhurst'                                                                                        
reference_long = 273                                                 
reference_lat = 0                                
plot_spirals = True                            
plot_sun_body_line = True                        
long_offset = 270                                
reference_vsw = 400                             
return_plot_object = False                      
transparent = False                            
numbered_markers = True                       

sm = SolarMACH(date, body_list, vsw_list, reference_long, reference_lat, coord_sys)

sm.plot(
   plot_spirals=plot_spirals,
   plot_sun_body_line=plot_sun_body_line,
   reference_vsw=reference_vsw,
   transparent=transparent,
   numbered_markers=numbered_markers,
   long_offset=long_offset,
   return_plot_object=return_plot_object,
   
)