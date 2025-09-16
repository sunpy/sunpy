.. _sunpy-how-to-mu:

****************************************************
Calculate the heliocentric angle and get :math:`\mu`
****************************************************

The heliocentric angle is the angle between direction from a point on the surface of the body toward the observer, with respect to the local vertical of that body.
The range of value is between 0 :math:`^\circ` (at the center) and 90 :math:`^\circ` (at the edge).

If one takes the :math:`\cos` of this value, this is referred to as :math:`\mu`.

It is worked out using the following equation:

.. math::

    heliocentric angle = \arcsin(\frac{d}{r})
    \mu=\cos(heliocentric angle)

where :math:`d` is the distance between the observer and the point on the solar surface and :math:`r` is the radius of the disk.

To calculate this using ``sunpy``, the function :func:`sunpy.coordinates.utils.get_heliocentric_angle` can be used.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> import numpy as np

    >>> from sunpy.coordinates.utils import get_heliocentric_angle

    >>> # At the disc center
    >>> hpc_coord_center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    >>> print(get_heliocentric_angle(hpc_coord_center))
    0.0 deg
    >>> # Mu
    >>> print(np.cos(get_heliocentric_angle(hpc_coord_center).to(u.radian)).to_value())
    1.0

    >>> # Almost at the limb
    >>> hpc_coord_limb = SkyCoord(944.35*u.arcsec, 0*u.arcsec, frame='helioprojective', observer="earth", obstime="2017-07-26")
    >>> print(get_heliocentric_angle(hpc_coord_limb))
    89.26429919 deg
    >>> # Mu
    >>> print(np.cos(get_heliocentric_angle(hpc_coord_limb).to(u.radian)).to_value())
    0.012840048535320435
