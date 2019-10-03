---
title: 'SunPy: A Python package for Solar Physics'
tags:
  - Python
  - Astronomy
  - Solar physics
authors:
  - name: Stuart J. Mumford
    orcid: 0000-0003-4217-4642
    affiliation: "1, 2"

affiliations:
 - name: SP2RC, School of Mathematics and Statistics, The University of Sheffield, UK
   index: 1
 - name: Aperio Software, Leeds, LS6 3HN, UK
   index: 2
date: 3 October 2019
bibliography: paper.bib
---

# Summary

The Sun, our nearest star, is a local laboratory for studying universal physical processes. Solar physics as a discipline includes studying the Sun both as a star and as the primary driver of space weather throughout the heliosphere. Due to the Sun's proximity, the temporal and spatial resolution of solar observations are orders of magnitude larger than those of other stars. This leads to significant differences in the data-analysis software needs of solar physicists compared with astrophysicists.

The `sunpy` Python package is a community-developed, free, and open-source solar data analysis environment for Python that provides a comprehensive set of tools for performing common tasks in solar data analysis. It is managed by the SunPy Project, an organization that facilitates and promotes the use of open development and open source packages like `sunpy` through community engagement and tools such as [GitHub](https://github.com/sunpy/sunpy/), [mailing lists](https://groups.google.com/forum/#!forum/sunpy), and [matrix](https://matrix.to/#/#sunpy:openastronomy.org).

The four most significant subpackages of `sunpy` are described below.

The `sunpy.net` subpackage provides a unified interface that simplifies and homogenizes search and retrieval by querying and downloading data from many solar data sources, irrespective of the underlying data-source client. It currently supports sourcing data from 18 different space- and ground-based solar observatories.

The `sunpy.map` and `sunpy.timeseries` subpackages provide core data types (`Map` and `TimeSeries`, respectively) that are designed to provide a general, standard, and consistent interface for loading and representing solar data across different instruments and missions. These classes load data which conform to solar physics standards and conventions such as FITS [@wells_fits_1981], FITS World Coordinate Systems (WCS) [@fits_wcs], and solar-specific FITS headers [@thompson_coordinates_2006], while allowing customization to account for differences in specific instruments. Visualization methods are also provided to inspect and plot those data. Example visulizations of both `TimeSeries` and `Map` are shown in Figure 1.

![map_timeseries_example](https://codimd.s3.shivering-isles.com/demo/uploads/upload_5eabd346e54c29b1d9b74aa11e351b30.png) *Left: An example of `TimeSeries` for the GOES X-ray Sensor in two broadband channels. Right: A `Map` of the extreme ultraviolet 171 $\mathring A$ channel of AIA corresponding to the time of a solar flare depicted by the vertical dashed line in the left-hand panel.*


The `sunpy.coordinates` subpackage provides support for representing and transforming coordinates used in solar physics and astrophysics. These coordinates may represent events (e.g., flares), features on or above the Sun (e.g., magnetic loops), or the position of structures traveling throughout the heliosphere (e.g., coronal mass ejections). The package currently implements the most widely used Sun-centered coordinate frames, and extends `astropy.coordinates`.

Other functionality provided by `sunpy` includes physical models of solar behavior, such as differential rotation, color maps for certain data sources, image-processing routines integrated with `Map`, and useful physical parameters such as constants.

The `sunpy` package is designed to be extensible, which means that it is easy to add support for additional instruments or data sources. It relies heavily on the `astropy` Python package as well as the scientific python stack (e.g. `numpy`, `scipy`, `matplotlib` and `pandas`).

A more complete description of the SunPy Project and the `sunpy` package, the methodology, development model, and implementation can be found in [@sunpy_community2019].

The SunPy Project supports affiliated packages, which build upon or extends the functionality of `sunpy`. Current affiliated packages are `drms` [@Glogowski2019drms], [`ndcube`](https://docs.sunpy.org/projects/ndcube), [`radiospectra`](https://docs.sunpy.org/projects/radiospectra) and [`IRISPy`](https://docs.sunpy.org/projects/irispy). The Project is also a member of the Python in Heliophysics community [PyHC, @annex], whose mission is to enable interdisciplinary analysis across all sub-disciplines of heliophysics by adhering to standards for code development and interoperability.


# Acknowledgements

SunPy is a NumFOCUS sponsored package.

# References
