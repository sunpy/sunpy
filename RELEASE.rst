The SunPy project is happy to announce the release of SunPy 1.0.0!

This release has been in development for 14 months, it is a large change from the 0.x series of SunPy releases with a lot of old functionality removed and some exciting new features.

The major highlights of this release are:

  - A complete transition of the whole code base to use `~astropy.time.Time`, which was implemented by Vishnunarayan K I as part of Google Summer of Code 2018.
  - A rewrite of how all the clients in `sunpy.net` download files from the internet. This means vastly improved progress bars, skipping downloads if files are present, and better visibility and retrying of download errors.
  - A rewrite of the differential rotation and image warping code to correctly account for observer location using the Astropy coordinate functionality.
  - Removal of many deprecated functions and submodules, we have used the 1.0 release as a chance to clean out SunPy reducing the number of lines of Python code in the project by almost 3,000!
  - The first release of SunPy to be Python 3 only, requiring Python 3.6+.

See `What's New in SunPy 1.0 <https://docs.sunpy.org/en/stable/whatsnew/1.0.html>`__ for more details and the `Full Changelog <https://docs.sunpy.org/en/stable/whatsnew/changelog.html>`__ for the full list of over 150 changes in SunPy 1.0.


This release includes many breaking changes and may require your code to be updated to support it.
We hope you see the value in having the extra features these changes enabled and a code base that is easier to maintain.
We will be continuing to release bug fixes to the 0.9.z series until the end of 2019.
We hope this gives you enough time to update your code and start enjoying all the improvements in SunPy 1.0.


By the numbers this release of SunPy contains 1913 commits in 332 merged pull requests closing 582 issues from 46 people, 25 of which are first time contributors to SunPy.

The people who have contributed to the code for this release are:

    Abhigyan Bose  *
    Airmansmith97  *
    Albert Y. Shih
    Andre Chicrala  *
    Andrew Hill  *
    Andrew Inglis
    Andrew Leonard
    Arthur Eigenbrot  *
    Bernhard M. Wiedemann  *
    Brandon Stone  *
    Brigitta Sipocz
    Daniel D'Avella
    Daniel Ryan
    David Pérez-Suárez
    David Stansby
    Deepankar Sharma  *
    Emmanuel Arias  *
    Harsh Mathur  *
    Jack Ireland
    Jacob  *
    Jai Ram Rideout  *
    Larry Manley
    Laura Hayes  *
    Manas Mangaonkar  *
    Matthew Mendero  *
    Md Akramul Haque  *
    Michael S Kirk
    Mickaël Schoentgen  *
    Monica Bobra  *
    Nabil Freij
    Naman9639  *
    Prateek Chanda
    Punyaslok Pattnaik
    Reid Gomillion  *
    Sarthak Jain  *
    Shane Maloney
    Shresth Verma
    Sourav Ghosh  *
    Steven Christe
    Stuart Mumford
    Vishnunarayan K I.
    Will Barnes
    Yash Jain
    Yash Krishan  *
    Yash Sharma  *
    jamescalixto  *

Where a * indicates their first contribution to SunPy.
