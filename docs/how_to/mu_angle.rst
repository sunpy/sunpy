.. _sunpy-how-to-mu-angle:

*************************************
Calculate the heliocentric (mu) angle
*************************************

`The following is based on a Stack Exchange answer <https://astronomy.stackexchange.com/a/48172>`__.

The heliocentric angle is the angle between direction from a point on the surface of the body toward the observer, with respect to the local vertical of that body.
The range of value is between 0 :math:`^\circ` (at the center) and 90 :math:`^\circ` (at the edge).
If one takes the :math:`\cos` of this value, this is referred to as :math:`\mu`.

It is worked out using the following equation:

.. math::

    heliocentric angle = \arcsin(\frac{d}{r})
    \mu=\cos(heliocentric angle)

where :math:`d` is the distance between the observer and the point on the solar surface and :math:`r` is the radius of the disk.

To calculate this using ``sunpy``, the functoion :func:`sunpy.coordinates.utils.get_mu_angle` can be used.

.. code-block:: python

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    from sunpy.coordinates.utils import get_mu_angle

    # At the disc center
    hpc_coord_center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    print(get_mu_angle(hpc_coord_center))

    # Almost at the limb
    hpc_coord_limb = SkyCoord(944.35*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    print(get_mu_angle(hpc_coord_limb))
