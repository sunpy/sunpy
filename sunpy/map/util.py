import sunpy
from sunpy import wcs
from scipy.interpolate import griddata
import numpy as np
#from sim.wave2d.wave2d import euler_zyz
#from matplotlib import colors

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

def map_hpc_to_hg(map, xbin = 1, ybin = 1):
	"""Take a map (like an AIA map) and convert it from Helio-projective Cartesian coordinates to Heliographic coordinates."""

	width = np.array(map.shape[1])
	height = np.array(map.shape[0])
	scale = np.array([map.scale.get('x'), map.scale.get('y')])
	crpix = np.array([map.reference_pixel.get('x'), map.reference_pixel.get('y')])
	crval = np.array([map.reference_coordinate.get('x'), map.reference_coordinate.get('y')])
	coordinate_system = [map.coordinate_system.get('x'), map.coordinate_system.get('y')]
	x,y = wcs.convert_pixel_to_data(width, height, scale[0], scale[1], crpix[0], crpix[1], crval[0], crval[1], coordinate_system[0])
	
	rsun = map.rsun_meters
	dsun = map.dsun
	
	b0 = map.heliographic_latitude
	l0 = map.heliographic_longitude
	units = [map.units.get('x'), map.units.get('y')]
	
	lon_map, lat_map = wcs.convert_hpc_hg(rsun, dsun, units[0], units[1], b0, l0, x, y)
	
	lon_bin = xbin
	lat_bin = ybin 
	lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
	lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))
	
	lon = np.arange(lon_range[0], lon_range[1], lon_bin)
	lat = np.arange(lat_range[0], lat_range[1], lat_bin)
	newgrid = np.meshgrid(lon, lat)
	
	points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
	values = np.array(map).ravel()
	
	# get rid of all of the bad (nan) indices (i.e. those off of the sun)
	index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
	points = np.vstack((points[index,0], points[index,1])).T
	
	values = values[index]
		
	newdata = griddata(points, values, newgrid, method="linear")
	
	header = map.get_header()
	header['CDELT1'] = lon_bin
	header['NAXIS1'] = len(lon)
	header['CRVAL1'] = lon.min()
	header['CRPIX1'] = 1
	header['CRPIX2'] = 1
	header['CUNIT1'] = "deg"
	header['CTYPE1'] = "HG"
	header['CDELT2'] = lat_bin
	header['NAXIS2'] = len(lat)
	header['CRVAL2'] = lat.min()
	header['CUNIT2'] = "deg"
	header['CTYPE2'] = "HG"
	
	transformed_map = sunpy.map.make_map(newdata, header)
	
	transformed_map.cmap = map.cmap
	transformed_map.name = map.name
	transformed_map.date = map.date
	
	return transformed_map

def map_hg_to_hpc(map, xbin = 10, ybin = 10):
    """Take a map in heliographic coordinates (HG) and convert it to 
    helioprojective cartesian coordinates (HPC)."""

    lon,lat = wcs.convert_pixel_to_data(map.header)
    x_map, y_map = wcs.convert_hg_hpc(map.header, lon, lat, units ='arcsec')
    
    x_range = (np.nanmin(x_map), np.nanmax(x_map))
    y_range = (np.nanmin(y_map), np.nanmax(y_map))
    
    x = np.arange(x_range[0], x_range[1], xbin)
    y = np.arange(y_range[0], y_range[1], ybin)
    newgrid = np.meshgrid(x, y)
    
    # newgrid = wcs.convert_hpc_hg(map.header, xgrid/(3600), ygrid/(3600))
    
    points = np.vstack((x_map.ravel(), y_map.ravel())).T
    values = np.array(map).ravel()
    newdata = griddata(points, values, newgrid, method="linear")
    
    # now grab the original map header and update it
    header = map.header.copy()
    header["CDELT1"] = xbin
    header["NAXIS1"] = len(x)
    header["CRVAL1"] = x.min()
    header["CRPIX1"] = 1
    header["CUNIT1"] = "arcsec"
    header["CTYPE1"] = "HPLT-TAN"
    header["CDELT2"] = ybin
    header["NAXIS2"] = len(y)
    header["CRVAL2"] = y.min()
    header["CRPIX2"] = 1
    header["CUNIT2"] = "arcsec"
    header["CTYPE2"] = "HPLT-TAN"
    
    transformed_map = sunpy.map.BaseMap(newdata, header)
    
    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date

    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}

    return transformed_map

def map_hpc_to_hg_rotate(map, epi_lon = 0, epi_lat = 0, xbin = 1, ybin = 1):
    """Take a map (like an AIA map) and convert it from HPC to HG."""

    #import sunpy
    #import util
    #from sunpy import wcs
    #import numpy as np
    #from scipy.interpolate import griddata
    from sim.wave2d.wave2d import euler_zyz
    #from matplotlib import colors
    
    # epi_lon = -10
    # epi_lat = 0
    
    #aia = sunpy.Map(sunpy.AIA_171_IMAGE).resample([500,500])
    # tmap = util.map_hpc_to_hg(aia)
    # tmap.show()
    
    #map = aia
    
    x, y = wcs.convert_pixel_to_data(map.header)
    
    hccx, hccy, hccz = wcs.convert_hpc_hcc_xyz(map.header, x, y)
    
    # rot_hccz, rot_hccy, rot_hccx = euler_zyz((hccz, hccx, hccy), (epi_lon, 90.-epi_lat, 0.))
    rot_hccz, rot_hccx, rot_hccy = euler_zyz((hccz, hccx, hccy), (0., epi_lat-90., -epi_lon))
    # zpp, xpp, ypp = euler_zyz(zxy_p, (0., hglt_obs, total_seconds*rotation))

    lon_map, lat_map = wcs.convert_hcc_hg(map.header, rot_hccx, rot_hccy, z = rot_hccz)
    
    lon_bin = xbin
    lat_bin = ybin 
    lon_range = (np.nanmin(lon_map), np.nanmax(lon_map))
    lat_range = (np.nanmin(lat_map), np.nanmax(lat_map))
    
    lon = np.arange(lon_range[0], lon_range[1], lon_bin)
    lat = np.arange(lat_range[0], lat_range[1], lat_bin)
    newgrid = np.meshgrid(lon, lat)
    
    #This extra conversion and rotation back are needed to determine where to
    #mask out points that can't have corresponding data
    ng_xyz = wcs.convert_hg_hcc_xyz(map.header, newgrid[0], newgrid[1])
    ng_zp, ng_xp, ng_yp = euler_zyz((ng_xyz[2], ng_xyz[0], ng_xyz[1]),
                                    (epi_lon, 90.-epi_lat, 0.))
    
    
    points = np.vstack((lon_map.ravel(), lat_map.ravel())).T
    values = np.array(map).ravel()
        
    # get rid of all of the bad (nan) indices (i.e. those off of the sun)
    index = np.isfinite(points[:,0]) * np.isfinite(points[:,1])
    #points = np.vstack((points[index,0], points[index,1])).T
    points = points[index]
    values = values[index]
    
    newdata = griddata(points, values, newgrid, method="cubic")
    newdata[ng_zp < 0] = np.nan
    
    header = map.header.copy()
    header['CDELT1'] = lon_bin
    header['NAXIS1'] = len(lon)
    header['CRVAL1'] = lon.min()
    header['CRPIX1'] = 1
    header['CRPIX2'] = 1
    header['CUNIT1'] = "deg"
    header['CTYPE1'] = "HG"
    header['CDELT2'] = lat_bin
    header['NAXIS2'] = len(lat)
    header['CRVAL2'] = lat.min()
    header['CUNIT2'] = "deg"
    header['CTYPE2'] = "HG"
    
    transformed_map = sunpy.map.BaseMap(newdata, header)
    
    transformed_map.cmap = map.cmap
    transformed_map.name = map.name
    transformed_map.date = map.date
    transformed_map.center = {
        "x": wcs.get_center(header, axis='x'),
        "y": wcs.get_center(header, axis='y')}
    
    #transformed_map.show(norm = colors.Normalize(0,1000))
    
    return transformed_map