import astropy.units as u
from astropy.coordinates import HeliocentricMeanEcliptic, SphericalRepresentation
import sunpy.coordinates.frames as f


class Transformation:
    """
    Handles transformations between different solar and astronomical coordinate frames.
    """

    frame_names = [
        'HGS', 'HGC', 'HCC', 'HPC', 'HCI', 
        'HAE', 'HEE', 'GSE', 'GEI'
    ]
    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    def setup_cache(self):
        """
        Initializes the cache of coordinate frames with predefined representations.

        Returns:
            dict: A dictionary mapping frame names to their respective frame objects.
        """
        obstime = '2023-01-01'
        vect = SphericalRepresentation(0 * u.deg, 0 * u.deg, 1 * u.AU)
        observer = f.HeliographicStonyhurst(vect, obstime=obstime)
        
        frames = {
            'HGS': f.HeliographicStonyhurst(vect, obstime=obstime),
            'HGC': f.HeliographicCarrington(vect, obstime=obstime, observer=observer),
            'HCC': f.Heliocentric(vect, obstime=obstime, observer=observer),
            'HPC': f.Helioprojective(vect, obstime=obstime, observer=observer),
            'HCI': f.HeliocentricInertial(vect, obstime=obstime),
            'HAE': HeliocentricMeanEcliptic(vect, obstime=obstime, equinox='J2000'),
            'HEE': f.HeliocentricEarthEcliptic(vect, obstime=obstime),
            'GSE': f.GeocentricSolarEcliptic(vect, obstime=obstime),
            'GEI': f.GeocentricEarthEquatorial(vect, obstime=obstime, equinox='J2000'),
        }
        return frames

    def setup(self, frames, src, dest):
        """
        Validates the source and destination frames for a transformation.

        Args:
            frames (dict): Dictionary of frames.
            src (str): Source frame name.
            dest (str): Destination frame name.

        Raises:
            NotImplementedError: If the source and destination frames are identical.
        """
        if src == dest:
            raise NotImplementedError("Source and destination frames must be different.")

    def time_transform(self, frames, src, dest):
        """
        Performs a transformation from the source frame to the destination frame.

        Args:
            frames (dict): Dictionary of frames.
            src (str): Source frame name.
            dest (str): Destination frame name.
        """
        frames[src].transform_to(frames[dest])
