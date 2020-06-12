The SunPy project is happy to announce the release of SunPy 2.0!
SunPy is an open-source Python library for Solar Physics data analysis and visualization.

This release is our second long(er) term support release, that we will be supporting with bug fixes until 3.0 in roughly a years time.
With this release the 1.0 and 1.1 releases will no longer recieve bug fixes and we encourage everyone to upgrade to 2.0

The major highlights of this release are:

* `~sunpy.net.Fido` now supports tab completion of search attributes.
  This allows you to do ``a.Instrument.AIA``, and print ``a.Instrument`` to see the list of known supported instruments.
* `~sunpy.instr.aia.aiaprep` has been deprecated in favour of the functionality in the `aiapy <https://aiapy.readthedocs.io/>`__ package.
* Various fixes and clarifications to pixel indexing in the map subpackage.
* Standidisation of specifying rectangles in coordinate space in the ``submap`` and ``draw_rectangle`` methods of ``GenericMap``.
* A HTML quicklook preview of ``GenericMap`` and ``MapSequence`` which can be accessed with the new ``quicklook()`` method.
  This is also the default display in Jupyter notebooks.
* Integration of differential rotation into the sunpy coordinate framework.
  This enables, amongst other things, the warping of images with the ``reproject`` package and the plotting of rotated grid lines with ``WCSAxes``.

See `What's New in SunPy 2.0 <https://docs.sunpy.org/en/stable/whatsnew/2.0.html>`__ for more details and the `Full Changelog <https://docs.sunpy.org/en/stable/whatsnew/changelog.html>`__ for the full list of over 100 changes in 2.0.


This release of SunPy contains 1044 commits in 290 merged pull requests closing 144 issues from 33 people, 16 of which are first time contributors to SunPy.

The people who have contributed to the code for this release are:

    Abhijeet Manhas  *
    Abijith B  *
    Albert Y. Shih
    Amogh J  *
    Arfon Smith  *
    Arib Alam  *
    David Pérez-Suárez
    David Stansby
    Deepankar Sharma
    Jack Ireland
    Jai Ram Rideout
    James Paul Mason
    Kris Akira Stern  *
    Laura Hayes
    Lazar Zivadinovic  *
    Mark Cheung  *
    Monica Bobra
    Nabil Freij
    Ole Streicher
    Pankaj Mishra  *
    Raahul Singh
    Rajiv Ranjan Singh
    Rutuja Surve  *
    Sarthak Jain
    Sashank Mishra  *
    Steven Christe
    Stuart Mumford
    Swapnil Kannojia  *
    Utkarsh Parkhi  *
    Will Barnes
    abijith-bahuleyan  *
    honey  *
    mridulpandey  *

Where a * indicates their first contribution to SunPy.
