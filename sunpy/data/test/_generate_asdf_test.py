"""
This file will generate an asdf file in the test data, using the newest schema version.
"""
import asdf
import astropy.units as u

from sunpy.data.test import rootdir


def generate_asdf_tree(obj, filename):
    tree = {"object": obj}
    with asdf.AsdfFile(tree) as af:
        # TODO: Automatically determine filename based on tag used
        af.write_to(rootdir / filename)


if __name__ == "__main__":
    import sunpy.map
    test_map = rootdir / "aia_171_level1.fits"
    obj = sunpy.map.Map(test_map)
    obj = obj.resample((2, 2)*u.pix)
    # obj = obj.shift(10*u.arcsec, 0*u.arcsec)
    generate_asdf_tree(obj, "aiamap_genericmap_1.0.0.asdf")
