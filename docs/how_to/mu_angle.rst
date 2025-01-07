.. _sunpy-how-to-mu-angle:

********************************
Calculate the heliocentric angle
********************************

`The following is based on a Stack Exchange answer <https://astronomy.stackexchange.com/a/48172
>__`.

The heliocentric angle is the angle between direction from a point on the surface of the body toward the observer, with respect to the local vertical of that body.
The range of value is between 0:math:`^\circ` (at the center) and 90:math:`^\circ` (at the edge).
If one takes the :math:`\cos` of this value, this is referred to as :math:`\mu`.

It is worked out using the following equation:

.. math::

    heliocentric angle = \arcsin(\frac{d}{r})
    \mu=\cos(heliocentric angle)

where :math:`d` is the distance between the observer and the point on the solar surface and :math:`r` is the radius of the disk.

To calculate this using ``sunpy``, the following code can be used:

.. code-block:: python

    from astropy.coordinates import SkyCoord
    from astropy.coordinates.representation import CartesianRepresentation
    import sunpy.coordinates
    import astropy.units as u
    import numpy as np

    # Disc center
    hpc_coord_center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    hcc = hpc_coord_center.heliocentric
    normal = hcc.cartesian
    to_observer = CartesianRepresentation(0, 0, 1) * hcc.observer.radius - normal
    mu = normal.dot(to_observer) / (normal.norm() * to_observer.norm())
    print(np.rad2deg(np.arccos(mu)))

    # Almost the limb
    hpc_coord_limb = SkyCoord(944.35*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    hcc = hpc_coord_limb.heliocentric
    normal = hcc.cartesian
    to_observer = CartesianRepresentation(0, 0, 1) * hcc.observer.radius - normal
    mu = normal.dot(to_observer) / (normal.norm() * to_observer.norm())
    print(np.rad2deg(np.arccos(mu)))
