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
  - name: Nabil Freij
  - name: Steven Christe
  - name: Jack Ireland
  - name: Florian Mayer
  - name: Keith Hughitt
  - name: Albert Y. Shih
  - name: Daniel Ryan
  - name: Simon Liedtke
  - name: David Pérez-Suárez
  - name: Pritish Chakraborty
  - name: Vishnunarayan K I.
  - name: Andrew Inglis
  - name: Punyaslok Pattnaik
  - name: Brigitta Sipocz
  - name: Rishabh Sharma
  - name: Andrew Leonard
  - name: David Stansby
  - name: Russell Hewett
  - name: Alex Hamilton
  - name: Laura Hayes
  - name: Asish Panda
  - name: Matt Earnshaw
  - name: Nitin Choudhary
  - name: Ankit Kumar
  - name: Prateek Chanda
  - name: Md Akramul Haque
  - name: Michael S Kirk
  - name: Michael Mueller
  - name: Sudarshan Konge
  - name: Rajul Srivastava
  - name: Yash Jain
  - name: Samuel Bennett
  - name: Ankit Baruah
  - name: Will Barnes
  - name: Michael Charlton
  - name: Cubostar
  - name: Shane Maloney
  - name: Nicky Chorley
  - name: Himanshu
  - name: Sanskar Modi
  - name: James Paul Mason
  - name: Naman9639
  - name: Yash
  - name: Jose Ivan Campos Rozo
  - name: Larry Manley
  - name: Agneet Chatterjee
  - name: John Evans
  - name: Michael Malocha
  - name: Monica Bobra
  - name: Sourav Ghosh
  - name: Airmansmith97
  - name: Dominik Stańczak
  - name: Ruben De Visscher
  - name: Shresth Verma
  - name: Ankit Agrawal
  - name: Dumindu Buddhika
  - name: Swapnil Sharma
  - name: Jongyeob Park
  - name: Matt Bates
  - name: Dhruv Goel
  - name: Garrison Taylor
  - name: Goran Cetusic
  - name: Jacob
  - name: Mateo Inchaurrandieta
  - name: Sally Dacie
  - name: Sanjeev Dubey
  - name: Deepankar Sharma
  - name: Erik M. Bray
  - name: Jai Ram Rideout
  - name: Serge Zahniy
  - name: Tomas Meszaros
  - name: Abhigyan Bose
  - name: Andre Chicrala
  - name: Ankit
  - name: Chloé Guennou
  - name: Daniel D'Avella
  - name: Daniel Williams
  - name: Jordan Ballew
  - name: Nick Murphy
  - name: Priyank Lodha
  - name: Thomas Robitaille
  - name: Yash Krishan
  - name: Andrew Hill
  - name: Arthur Eigenbrot
  - name: Benjamin Mampaey
  - name: Bernhard M. Wiedemann
  - name: Carlos Molina
  - name: Duygu Keşkek
  - name: Ishtyaq Habib
  - name: Joe Letts
  - name: Juanjo Bazán
  - name: Quinn Arbolante
  - name: Reid Gomillion
  - name: Yash Kothari
  - name: Yash Sharma
  - name: Abigail Stevens
  - name: Adrian Price-Whelan
  - name: Ambar Mehrotra
  - name: Arseniy Kustov
  - name: Brandon Stone
  - name: Dang Trung Kien
  - name: Emmanuel Arias
  - name: Fionnlagh Mackenzie Dover
  - name: Freek Verstringe
  - name: Gulshan Mittal
  - name: Harsh Mathur
  - name: Igor Babuschkin
  - name: Jaylen Wimbish
  - name: Juan Camilo Buitrago-Casas
  - name: Kalpesh Krishna
  - name: Kaustubh Hiware
  - name: Manas Mangaonkar
  - name: Matthew Mendero
  - name: Mickaël Schoentgen
  - name: Norbert Gyenge
  - name: Ole Streicher
  - name: Rajasekhar Reddy Mekala
  - name: Rishabh Mishra
  - name: S Shashank
  - name: Sarthak Jain
  - name: Tannmay Yadav
  - name: Tessa D. Wilkinson
  - name: Tiago Pereira
  - name: Yudhik Agrawal
  - name: jamescalixto
  - name: yasintoda

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

![Left: An example of `TimeSeries` for the GOES X-ray Sensor in two broadband channels. Right: A `Map` of the extreme ultraviolet 171 $\mathring A$ channel of AIA corresponding to the time of a solar flare depicted by the vertical dashed line in the left-hand panel.](https://codimd.s3.shivering-isles.com/demo/uploads/upload_5eabd346e54c29b1d9b74aa11e351b30.png)

The `sunpy.coordinates` subpackage provides support for representing and transforming coordinates used in solar physics and astrophysics. These coordinates may represent events (e.g., flares), features on or above the Sun (e.g., magnetic loops), or the position of structures traveling throughout the heliosphere (e.g., coronal mass ejections). The package currently implements the most widely used Sun-centered coordinate frames, and extends `astropy.coordinates`.

Other functionality provided by `sunpy` includes physical models of solar behavior, such as differential rotation, color maps for certain data sources, image-processing routines integrated with `Map`, and useful physical parameters such as constants.

The `sunpy` package is designed to be extensible, which means that it is easy to add support for additional instruments or data sources. It relies heavily on the `astropy` Python package as well as the scientific python stack (e.g. `numpy`, `scipy`, `matplotlib` and `pandas`).

A more complete description of the SunPy Project and the `sunpy` package, the methodology, development model, and implementation can be found in [@sunpy_community2019].

The SunPy Project supports affiliated packages, which build upon or extends the functionality of `sunpy`. Current affiliated packages are `drms` [@Glogowski2019drms], [`ndcube`](https://docs.sunpy.org/projects/ndcube), [`radiospectra`](https://docs.sunpy.org/projects/radiospectra) and [`IRISPy`](https://docs.sunpy.org/projects/irispy). The Project is also a member of the Python in Heliophysics community [PyHC, @annex], whose mission is to enable interdisciplinary analysis across all sub-disciplines of heliophysics by adhering to standards for code development and interoperability.


# Acknowledgements

SunPy is a NumFOCUS sponsored package.

# References
