from sunpy.extern import six
import sunpy.map

def sswmap2sunpymap(sswmap):
    '''
    This is a helpper function to convert solarsoft maps restored using readsav
    
    It's been created using the values as defined in `make_map.pro` from
    http://hesperia.gsfc.nasa.gov/ssw/gen/idl/maps/make_map.pro
    On 24 Oct 2015 (latest modification 22 Sep 2014 by Zarro)
    '''

    # todo Cube?

    # --- IDL code ---
    # ap=create_struct('time',stime,'id',id,'dur',dur,$
    #                  'xunits',xunits,'yunits',yunits,$
    #                  'roll_angle',double((roll_angle mod 360.)),$
    #                  'roll_center',reform(double(roll_center)))

    # VMS format is the default of a not valid time: 18-Jan-1988 17:20:43.123
    vms_fmt = '%d-%b-%Y %H:%M:%S.%f'
    time =  datetime.datetime.strptime(sswmap['time'][0], vms_fmt)
    # This may not be the obs_date, but the date to which a map has been rotated

    # id; unique string identifier. TODO: Check index2map to see what it is
    id_map = sswmap['id'][0]

    # dur: TODO: Check, exp_time duration?
    dur = sswmap['dur'][0] # defaults to 0

    # x/yunits -> cunit1/2
    # If the units are not a string (or empty) defaults to arcsecs
    axis_units = {'cunit1': 'xunits', 'cunit2': 'yunits'}
    for key, value in six.iteritems(axis_units):
        axis_units[key] = sswmap[value][0]
        # arcsecs is not a valid recognised unit.
        if axis_units[key] == 'arcsecs' then axis_units[key] = 'arcsec' 
    # TODO: Check how we handle empty or other units used by plot_map

    # roll_angle: image roll (deg clockwise from N) - defaults to 0.
    roll_angle = sswmap['roll_angle'][0] # deg => TODO: check how our map interprets it

    # roll_center: optional and if not set it defaults to [xc, yc].
    roll_center = sswmap['roll_center'][0] # 2-element array

    
    # --- IDL code ---
    # if old_format then map=create_struct('xp',xp,'yp',yp,map) else $
    # map=create_struct('xc',xc,'yc',yc,'dx',dx,'dy',dy,map)

    center_coordinates = {'crval1': ['xc', 'xp'], 'crval2': ['yc', 'yp']}
    for key, value in six.iteritems(center_coordinates):
        center_coordinates[key] = sswmap.get(value[0], sswmap[value[1]])[0]
    # TODO: Check whether xp is the same that xc or other kind of array

    scale_pixels = {'cdelt1': 'dx', 'cdelt2': 'dy'}
    for key, value in six.iteritems(scale_pixels):
        scale_pixels[key] = sswmap[value][0]
        
    # --- IDL code ---
    # map=create_struct('data',data,map)
    data = sswmap['data'][0]

    # --- IDL code ---
    # if is_struct(extra) then map=create_struct(map,extra)

    # extra?; TODO: Check what can be in here
    # date_obs? ctype? r_sun? telescope?

    header = dict(center_coordinates, **scale_pixels}
    header.update(axis_units)
    header.update({'id': id_map, 'time': time, 'dur': dur,
                   'roll_angle': roll_angle, 'roll_center': roll_center})
    # TODO: Check whether these are the appropriate keywords for us
    
    # TODO: what are the minimum we need for a map?
    return sunpy.map.Map(data, header)
