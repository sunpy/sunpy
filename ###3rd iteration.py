###3rd iteration

import numpy as np
import rasterio

'''
I have used SLDEM2015 data for the DEM data. The data is available at the following link: https://pgda.gsfc.nasa.gov/products/54
I have only used their data to process the .j2 files for further usage
'''

def process_j2_dem(file_path):
    """
    I found .j2 DEM files from the official NASA website. hence considered using .j2 files for DEM data.
    Processing a .j2 DEM file to extract the latitude, longitude, and elevation data.

    Parameters
    ----------
    file_path : str
        Path to the .j2 DEM file.

    Returns
    -------
    dem_data : list of tuples
        A list of (latitude, longitude, elevation) tuples for each pixel in the DEM.
    """
    # Open the .j2 DEM file using rasterio
    with rasterio.open(file_path) as src:
        # Read the elevation data
        elevation_data = src.read(1)  # First band contains elevation values
        
        # Get the spatial resolution and bounds
        transform = src.transform
        width, height = src.width, src.height

        # Create arrays for latitude and longitude
        lon_min, lat_max = transform * (0, 0)  # Top-left corner
        lon_max, lat_min = transform * (width, height)  # Bottom-right corner

        # Create latitude and longitude arrays
        lons = np.linspace(lon_min, lon_max, width)
        lats = np.linspace(lat_max, lat_min, height)

        # Convert the 2D elevation array into a list of (lat, lon, elevation)
        dem_data = []
        for i, lat in enumerate(lats):
            for j, lon in enumerate(lons):
                elevation = elevation_data[i, j]
                # Only include valid elevation values (filter out no-data values)
                if not np.isnan(elevation):
                    dem_data.append((lat, lon, elevation))

    return dem_data


def compute_limb_profile_from_dem(dem_data, observer_distance, resolution=18000):
    """
    Compute the lunar limb profile from DEM data, with integrated transformations.

    Parameters
    ----------
    dem_data : array
        Lunar DEM data, formatted as [(latitude, longitude, elevation), ...].
    observer_distance : float
        Distance of the observer from the Moon in km.
    resolution : int, optional
        Number of angular intervals for the limb profile (default is 18000).

    Returns
    -------
    limb_profile : np.ndarray
        Array of angular radii (α) at evenly spaced azimuthal angles.
    """
    # Initialize the limb profile
    limb_profile = np.zeros(resolution)

    # Precompute azimuthal bins
    angular_intervals = np.linspace(0, 2 * np.pi, resolution, endpoint=False)

    for lat, lon, elevation in dem_data:
        # Convert (lat, lon, elevation) to Cartesian coordinates
        lat, lon = np.radians(lat), np.radians(lon)
        R = 1737.4 + elevation  # Moon's radius + elevation in km
        x = R * np.cos(lat) * np.cos(lon)
        y = R * np.cos(lat) * np.sin(lon)
        z = R * np.sin(lat)

        # Apply a rotation for the observer's perspective
        x_rot, y_rot, z_rot = x, y, z  # Identity rotation for now (replace with real rotation)

        # Convert to polar coordinates relative to observer
        r = np.sqrt(y_rot**2 + z_rot**2)
        theta = np.arctan2(z_rot, y_rot)
        alpha = np.arctan(r / (observer_distance - x_rot))

        # Update the limb profile at the appropriate azimuthal index
        index = int((theta % (2 * np.pi)) / (2 * np.pi) * resolution)
        limb_profile[index] = max(limb_profile[index], alpha)

    return limb_profile


def eclipse_amount_with_dynamic_limb(observer, dem_data, observer_distance, *, resolution=18000):
    """
    Return the percentage of the Sun that is eclipsed by the Moon using a dynamic lunar limb model.

    Parameters
    ----------
    observer : `~astropy.coordinates.SkyCoord`
        The observer location and observation time.
    dem_data : array
        Lunar DEM data, formatted as [(latitude, longitude, elevation), ...].
    observer_distance : float
        Distance of the observer from the Moon in km.
    resolution : int, optional
        Number of angular intervals for the limb profile (default is 18000).

    Returns
    -------
    fraction : `~astropy.units.Quantity`
        The eclipse fraction as a percentage.
    """
    # Compute the lunar limb profile
    limb_profile = compute_limb_profile_from_dem(dem_data, observer_distance, resolution)

    # Observer's position relative to the Moon
    moon = get_body_heliographic_stonyhurst('moon', observer.obstime, observer=observer, quiet=True)
    observer = observer.transform_to(moon)

    # Compute azimuthal angle for limb profile
    vec_moon = moon.cartesian - observer.cartesian
    theta = np.arctan2(vec_moon.z, vec_moon.y)

    # Retrieve dynamic angular radius from limb profile
    index = int((theta % (2 * np.pi)) / (2 * np.pi) * resolution)
    m = limb_profile[index]

    # Compute Sun's angular radius and angular separation
    s = np.arcsin(constants.radius / observer.radius).value
    d = np.arccos(vec_sun.dot(vec_moon) / (observer.radius * vec_moon.norm())).value

    # Precalculate eclipse fraction using the dynamic lunar radius
    cs, ss = np.cos(s), np.sin(s)
    cm, sm = np.cos(m), np.sin(m)
    cd, sd = np.cos(d), np.sin(d)
    area_s = 2 * np.pi * (1 - cs)
    area_m = 2 * np.pi * (1 - cm)
    area_int = 2 * (np.pi
                    - np.arccos((cd - cs * cm) / (ss * sm))
                    - np.arccos((cm - cd * cs) / (sd * ss)) * cs
                    - np.arccos((cs - cd * cm) / (sd * sm)) * cm)

    # Handle edge cases for zero, total, and annular eclipses
    area_int[d >= s + m] = 0  # zero eclipse
    area_int[m >= s + d] = area_s[m >= s + d]  # total eclipse
    area_int[s >= m + d] = area_m[s >= m + d]  # annular eclipse

    # Divide by Sun's area and return fraction
    fraction = area_int / area_s
    fraction = np.clip(fraction, 0, 1)

    return u.Quantity(fraction.reshape(observer.data.shape)).to(u.percent)
