from asv_runner.benchmarks.mark import SkipNotImplemented

import astropy.units as u
from astropy.coordinates import HCRS, ITRS, HeliocentricMeanEcliptic, SphericalRepresentation

import sunpy.coordinates.frames as f


class TransformationHeliographic:
    frame_names = ['HCRS', 'HGS', 'HGC', 'HCC', 'HPC', 'HPR', 'HCI']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    def setup_cache(self):
        obstime = '2023-01-01'
        vect = SphericalRepresentation(0*u.deg, 0*u.deg, 1*u.AU)
        observer = f.HeliographicStonyhurst(vect, obstime=obstime)
        frames = {
            'HCRS': HCRS(vect, obstime=obstime),
            'HGS': f.HeliographicStonyhurst(vect, obstime=obstime),
            'HGC': f.HeliographicCarrington(vect, obstime=obstime, observer=observer),
            'HCC': f.Heliocentric(vect, obstime=obstime, observer=observer),
            'HPC': f.Helioprojective(vect, obstime=obstime, observer=observer),
            'HPR': f.HelioprojectiveRadial(vect, obstime=obstime, observer=observer),
            'HCI': f.HeliocentricInertial(vect, obstime=obstime),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        frames[src].transform_to(frames[dest])


class TransformationEcliptic:
    frame_names = ['HAE', 'HEE', 'GSE', 'GEI']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    def setup_cache(self):
        obstime = '2023-01-01'
        vect = SphericalRepresentation(0*u.deg, 0*u.deg, 1*u.AU)
        frames = {
            'HAE': HeliocentricMeanEcliptic(vect, obstime=obstime, equinox='J2000'),
            'HEE': f.HeliocentricEarthEcliptic(vect, obstime=obstime),
            'GSE': f.GeocentricSolarEcliptic(vect, obstime=obstime),
            'GEI': f.GeocentricEarthEquatorial(vect, obstime=obstime, equinox='J2000'),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        frames[src].transform_to(frames[dest])


class TransformationMagnetic:
    frame_names = ['GEO', 'MAG', 'SM', 'GSM']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    def setup_cache(self):
        obstime = '2023-01-01'
        vect = SphericalRepresentation(0*u.deg, 0*u.deg, 1*u.AU)
        frames = {
            'GEO': ITRS(vect, obstime=obstime),
            'MAG': f.Geomagnetic(vect, obstime=obstime),
            'SM': f.SolarMagnetic(vect, obstime=obstime),
            'GSM': f.GeocentricSolarMagnetospheric(vect, obstime=obstime),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        frames[src].transform_to(frames[dest])
