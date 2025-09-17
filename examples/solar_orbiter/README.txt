Solar Orbiter Examples
=====================

This folder contains examples of how to work with data from the Solar Orbiter mission using SunPy and related packages.

About Solar Orbiter
====================

Solar Orbiter is a space mission of international collaboration between the European Space Agency (ESA) and NASA. The mission is designed to study the Sun and its influence on the heliosphere from a unique perspective. Unlike Earth-orbiting solar observatories, Solar Orbiter follows an elliptical orbit around the Sun that brings it as close as 0.28 AU to the Sun, and its orbital inclination will eventually allow it to observe the Sun's polar regions.

Key Instruments
===============

The examples in this folder focus on data from Solar Orbiter's remote sensing instruments:

* **EUI (Extreme Ultraviolet Imager)**: Images the solar corona and transition region

  - FSI (Full Sun Imager): Images the full solar disk in 174 Å and 304 Å
  - HRI_EUV: High-resolution images in 174 Å
  - HRI_LYA: High-resolution images in Lyman-alpha (1216 Å)

* **STIX (Spectrometer/Telescope for Imaging X-rays)**: X-ray imaging and spectroscopy

* **Metis**: Solar coronagraph for visible and UV observations

* **SoloHI (Solar Orbiter Heliospheric Imager)**: Wide-field heliospheric imager

Data Access
===========

Solar Orbiter data can be accessed through the Solar Orbiter Archive (SOAR) using the `sunpy-soar <https://pypi.org/project/sunpy-soar/>`__ package. This package integrates with SunPy's `~sunpy.net.Fido` interface to provide seamless data discovery and download.

Installation
------------

To run these examples, you'll need to install the ``sunpy-soar`` package::

    pip install sunpy-soar

Or using conda::

    conda install -c conda-forge sunpy-soar

Examples Overview
=================

The examples in this folder demonstrate:

1. **downloading_solo_eui_data.py**: Basic data search and download using sunpy-soar
2. **plotting_eui_images.py**: Visualization techniques for EUI images with proper coordinate handling
3. **advanced_solo_analysis.py**: Multi-instrument analysis and orbital perspective considerations

Unique Scientific Value
=======================

Solar Orbiter's orbital characteristics enable unique scientific observations:

* **Variable heliocentric distance**: Observations from 0.28 to 1.4 AU provide different spatial scales
* **Out-of-ecliptic perspective**: Eventually up to 33° inclination for polar observations
* **Co-rotation with solar features**: Enables tracking of solar structures over multiple rotations
* **Stereoscopic observations**: Combined with Earth-based observatories for 3D reconstruction
* **In-situ and remote sensing coordination**: Simultaneous particle and field measurements with imaging

Data Characteristics
====================

Solar Orbiter data has several important characteristics to consider:

**Observer Location**
Solar Orbiter's position varies significantly throughout its orbit. Always check the observer coordinates in the metadata:

* ``HGLN_OBS``: Heliographic longitude of observer
* ``HGLT_OBS``: Heliographic latitude of observer
* ``DSUN_OBS``: Distance from Sun to observer

**Coordinate Systems**
Due to the variable observer location, coordinate transformations between Solar Orbiter's perspective and Earth-based observations require careful handling of the observer position.

**Data Levels**
- **Level 0**: Raw, uncalibrated data
- **Level 1**: Calibrated data with instrument corrections applied
- **Level 2**: Science-ready data products (recommended for most analyses)
- **Level 3**: Higher-level data products and derived quantities

Tips for Analysis
=================

1. **Always check the observer position** when interpreting coordinate information
2. **Use appropriate coordinate transformations** when comparing with Earth-based observations
3. **Consider the viewing angle** when identifying solar features
4. **Take advantage of the unique perspective** for studies requiring off-limb observations
5. **Coordinate with other missions** for multi-viewpoint studies

Useful Resources
================

* `Solar Orbiter Mission Website <https://sci.esa.int/web/solar-orbiter/>`__
* `SOAR Data Archive <https://soar.esac.esa.int/soar/>`__
* `SunPy Documentation <https://docs.sunpy.org/>`__
* `sunpy-soar Documentation <https://docs.sunpy.org/projects/sunpy-soar/>`__
* `Solar Orbiter Science Archive Manual <https://www.cosmos.esa.int/web/soar/science-archive-manual>`__

For questions about Solar Orbiter data or these examples, please refer to the SunPy community resources or the Solar Orbiter Science Operations Centre documentation.