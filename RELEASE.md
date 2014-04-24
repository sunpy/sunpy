The SunPy development team is pleased to announce the release of SunPy 0.4.0.
This release contains many new features, some of which are contributed by people 
who participated in GSOC 2013. It includes the addition of a new local database for storing
and searching data, it features a HEK to VSO translator and a new HELIO module in net.
As well as this major work has been undertaken on the documentation and a new website developed.

This release contains 1025 commits from 19 people.

New Features:

    * **Major** documentation refactor. A far reaching re-write and restructure.
    * Add a SunPy Database to store and search local data.
    * Add beta support for querying the HELIO HEC
    * Add beta HEK to VSO query translation.
    * Add the ability to download the GOES event list.
    * Add support for downloading and querying the LYTAF database.
    * Add support for ANA data.
    * Updated sun.constants to use astropy.constants objects which include units, source,
    and error instide. For more info check out http://docs.astropy.org/en/latest/constants/index.html
    * Add some beta support for IRIS data products
    * Add a new MapCubeAnimator class with interactive widgets which is returned by mapcube.peek().
    * The Glymur library is now used to read JPEG2000 files.
    * GOESLightCurve now supports all GOES satellites.



The people who have contributed to this release are:

    Stuart Mumford
    Simon Liedtke
    Steven Christe
    Jack Ireland
    Andrew Inglis
    Nabil Freij
    Samuel Bennett
    David Perez-Suarez
    Pritish Chakraborty
    Albert Y. Shih
    John Evans
    Michael Malocha
    Florian Mayer
    Russell Hewett
    Jose Iván Campos
    Keith Hughitt
    Tiago Pereira
