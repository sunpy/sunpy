import astropy.units as u
from astropy.coordinates import HeliocentricMeanEcliptic, SphericalRepresentation

import sunpy.coordinates.frames as f


class Transformation:
    frame_names = ['HGS', 'HGC', 'HCC', 'HPC', 'HCI', 'HAE', 'HEE', 'GSE', 'GEI']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    def setup_cache(self):
        obstime = '2023-01-01'
        vect = SphericalRepresentation(0*u.deg, 0*u.deg, 1*u.AU)
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
        if src == dest:
            raise NotImplementedError

    def time_transform(self, frames, src, dest):
        frames[src].transform_to(frames[dest])
