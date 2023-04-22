.. _how_to_search_the_vso

Searching the Virtual Solar Observatory (VSO)
=============================================

To search the Virtual Solar Observatory (VSO) for SDO AIA data in all channels over a given time range,
use the timerange (`~sunpy.net.attrs.Time`) and the instrument (`~sunpy.net.attrs.Instrument`) attrs,

.. code-block:: python

    >>> import astropy.units as u
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> time_range = a.Time('2020/03/04 00:00:15', '2020/03/04 00:00:30')
    >>> result = Fido.search(time_range, a.Instrument.aia)
    >>> result

You can also search for specific physical observables using `~sunpy.net.attr.Physobs`.
For example, you can search of line-of-sight (LOS) magnetic field measurements from HMI in
the same time range,

.. code-block:: python

    >>> Fido.search(time_range, a.Instrument.hmi, a.Physobs.los_magnetic_field)

You can also use relational operators when constructing queries.
For example, the AIA query above can also be expressed using the AND (&) operator,

.. code-block:: python

    >>> Fido.search(time_range & a.Instrument.aia)

Additionally, multiple operators can be chained together to, for example, search for only the 171 :math:`\AA` channel,

.. code-block:: python

    >>> Fido.search(time_range & a.Instrument.aia & a.Wavelength(171*u.angstrom))

The OR operator (|) can also be used to construct queries.
For example, to search for AIA data in this same time range from both the 94 and 171 channels,

.. code-block:: python

    >>> Fido.search(time_range,
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))

These relational operators are particularly useful when searching a given time interval for multiple instruments.
For example, to find the HMI LOS magnetic field data and the AIA 94 and 171 data in the given time interval,

.. code-block:: python

    >>> aia_params = a.Instrument.aia & (a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))
    >>> hmi_params = a.Instrument.hmi & a.Physobs.los_magnetic_field
    >>> Fido.search(time_range, aia_params | hmi_params)
