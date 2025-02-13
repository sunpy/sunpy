import numpy as np
from asv_runner.benchmarks.mark import SkipNotImplemented

import astropy.units as u
from astropy.coordinates import HCRS, ITRS, HeliocentricMeanEcliptic, SphericalRepresentation

import sunpy.coordinates.frames as f


class TransformationHeliographic:
    frame_names = ['HCRS', 'HGS', 'HGC', 'HCC', 'HPC', 'HPR', 'HCI']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    vect = SphericalRepresentation(np.arange(100001)*u.deg,
                                   np.linspace(-90, 90, 100001)*u.deg,
                                   np.linspace(0, 2, 100001)*u.AU)

    def setup_cache(self):
        obstime = '2023-01-01'
        observer = f.HeliographicStonyhurst(SphericalRepresentation(10*u.deg, 20*u.deg, 1*u.AU), obstime=obstime)
        frames = {
            'HCRS': HCRS(obstime=obstime),
            'HGS': f.HeliographicStonyhurst(obstime=obstime),
            'HGC': f.HeliographicCarrington(obstime=obstime, observer=observer),
            'HCC': f.Heliocentric(obstime=obstime, observer=observer),
            'HPC': f.Helioprojective(obstime=obstime, observer=observer),
            'HPR': f.HelioprojectiveRadial(obstime=obstime, observer=observer),
            'HCI': f.HeliocentricInertial(obstime=obstime),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        coord = frames[src].realize_frame(self.vect, copy=True)
        coord.transform_to(frames[dest])


class TransformationEcliptic:
    frame_names = ['HAE', 'HEE', 'GSE', 'GEI']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    vect = SphericalRepresentation(np.arange(100001)*u.deg,
                                   np.linspace(-90, 90, 100001)*u.deg,
                                   np.linspace(0, 2, 100001)*u.AU)

    def setup_cache(self):
        obstime = '2023-01-01'
        frames = {
            'HAE': HeliocentricMeanEcliptic(obstime=obstime, equinox='J2000'),
            'HEE': f.HeliocentricEarthEcliptic(obstime=obstime),
            'GSE': f.GeocentricSolarEcliptic(obstime=obstime),
            'GEI': f.GeocentricEarthEquatorial(obstime=obstime, equinox='J2000'),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        coord = frames[src].realize_frame(self.vect, copy=True)
        coord.transform_to(frames[dest])


class TransformationMagnetic:
    frame_names = ['GEO', 'MAG', 'SM', 'GSM']

    params = (frame_names, frame_names)
    param_names = ['src', 'dest']

    vect = SphericalRepresentation(np.arange(100001)*u.deg,
                                   np.linspace(-90, 90, 100001)*u.deg,
                                   np.linspace(0, 2, 100001)*u.AU)

    def setup_cache(self):
        obstime = '2023-01-01'
        frames = {
            'GEO': ITRS(obstime=obstime),
            'MAG': f.Geomagnetic(obstime=obstime),
            'SM': f.SolarMagnetic(obstime=obstime),
            'GSM': f.GeocentricSolarMagnetospheric(obstime=obstime),
        }
        return frames

    def setup(self, frames, src, dest):
        if src == dest:
            raise SkipNotImplemented

    def time_transform(self, frames, src, dest):
        coord = frames[src].realize_frame(self.vect, copy=True)
        coord.transform_to(frames[dest])
