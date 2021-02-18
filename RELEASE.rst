The SunPy project is happy to announce the release of SunPy 2.1!
SunPy is an open-source Python library for Solar Physics data analysis and visualization.

This is an intermediate release between the 2.0 and 3.0 long term support releases.

The major highlights of this release are:

* Support for the GOES-16 and GOES-17 XRS 1-second data. This data can now be queried and downloaded with `sunpy.net.Fido` and the files read in and analysed as a `sunpy.timeseries.TimeSeries`.
* ``sunpy.net.Fido`` now supports metadata only searches; currently supported metadata providers are HEK, HELIO and JSOC.
* ``sunpy.net.Fido`` now returns results to users in Astropy Table objects to allow easier manipulation and filtering of results.
* Cutouts of regions of interest can now be requested from the JSOC export service through ``sunpy.net.Fido``.
* It is now supported to transform coordinates with attached velocities, and the various ephemeris functions can optionally include velocity information. Transformations between coordinate frames will account for both any change in orientation of the velocity vector and any induced velocity due to relative motion between the frames.
* Several functions in `sunpy.map` have been significantly sped up with improved algorithms.


See `What's New in SunPy 2.1 <https://docs.sunpy.org/en/stable/whatsnew/2.1.html>`__ for more details and the `Full Changelog <https://docs.sunpy.org/en/stable/whatsnew/changelog.html>`__ for the full list of over 160 changes in 2.1.

This release of SunPy contains 1148 commits in 306 merged pull requests closing 202 issues from 38 people, 19 of which are first-time contributors to SunPy.

The people who have contributed to the code for this release are:

-  Abhijeet Manhas
-  Abhishek Pandey  *
-  Adrian Price-Whelan
-  Albert Y. Shih
-  Aryan Chouhan  *
-  Conor MacBride  *
-  Daniel Ryan
-  David Pérez-Suárez
-  David Stansby
-  Dipanshu Verma  *
-  Erik Tollerud  *
-  Jai Ram Rideout
-  Jeffrey Aaron Paul  *
-  Johan L. Freiherr von Forstner  *
-  Kateryna Ivashkiv  *
-  Koustav Ghosh  *
-  Kris Akira Stern
-  Kritika Ranjan  *
-  Laura Hayes
-  Lazar Zivadinovic
-  Nabil Freij
-  Rutuja Surve
-  Sashank Mishra
-  Shane Maloney
-  Shubham Jain  *
-  SophieLemos  *
-  Steven Christe
-  Stuart Mumford
-  Sudeep Sidhu  *
-  Tathagata Paul  *
-  Thomas A Caswell  *
-  Will Barnes
-  honey
-  mridulpandey
-  nakul-shahdadpuri  *
-  platipo  *
-  resakra  *
-  sophielemos  *

Where a * indicates that this release contains their first contribution to SunPy.
