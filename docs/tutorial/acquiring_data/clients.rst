.. doctest-skip-all

.. _sunpy-tutorial-acquiring-data-clients:

****************************
Overview of Supported Clients
****************************

``sunpy`` provides access to numerous remote solar physics data archives and web services through the `sunpy.net` sub-package.
While most users interact with `~sunpy.net.Fido`—the unified search and retrieval interface—queries are executed under the hood by specific **clients**.
Each client knows how to communicate with a particular archive, observatory, or data provider.

This page provides a comprehensive reference of all clients supported by `sunpy.net`, the missions and data products they offer, and the search attributes they require.

.. contents:: Contents
   :local:
   :depth: 2

.. _sunpy-tutorial-clients-fido-routing:

How Fido Routes Queries to Clients
==================================

When you call `~sunpy.net.Fido.search`, ``Fido`` inspects the search attributes you passed (such as `~sunpy.net.attrs.Time`, `~sunpy.net.attrs.Instrument`, or client-specific attributes) and automatically matches them to the appropriate client(s).

To see all registered clients available in your current installation, you can print ``Fido``:

.. code-block:: python

    >>> from sunpy.net import Fido
    >>> print(Fido)

You can also inspect registered client classes programmatically:

.. code-block:: python

    >>> from sunpy.net.base_client import BaseClient
    >>> for client in BaseClient._registry.keys():
    ...     print(client.__name__)

If a query matches multiple clients (or if you want to target a specific data source), you can restrict your search by supplying client-specific attributes or using the specific client class directly.

Summary of Supported Clients
============================

The following table summarizes the data clients available in ``sunpy.net``:

.. list-table:: Supported Data Clients in ``sunpy.net``
   :header-rows: 1
   :widths: 18 22 40 20

   * - Client Name
     - Module Path
     - Datasets & Instruments Provided
     - Key Search Attributes
   * - **VSOClient**
     - `~sunpy.net.vso.VSOClient`
     - Broad repository for solar images and spectral data (SDO, SOHO, STEREO, Hinode, BBSO, Kanzelhöhe, etc.).
     - `a.Time`, `a.Instrument`, `a.Source`, `a.Provider`, `a.Wavelength`, `a.Physobs`
   * - **JSOCClient**
     - `~sunpy.net.jsoc.JSOCClient`
     - Full-cadence SDO (AIA, HMI) data series, cutouts, Dopplergrams, and magnetograms from Stanford JSOC.
     - `a.Time`, `a.jsoc.Series`, `a.jsoc.Segment`, `a.jsoc.Notify`
   * - **SOARClient**
     - `~sunpy.net.soar.SOARClient`
     - All remote sensing (EUI, PHI, SPICE, STIX, Metis, SoloHI) and in-situ instruments on ESA/NASA Solar Orbiter.
     - `a.Time`, `a.Instrument`, `a.Level`, `a.soar.Product`
   * - **CDAWEBClient**
     - `~sunpy.net.cdaweb.CDAWEBClient`
     - In-situ solar wind, magnetic field, and particle data from NASA CDAWeb (Parker Solar Probe, Wind, ACE, etc.).
     - `a.Time`, `a.cdaweb.Dataset`
   * - **SOLARNETClient**
     - `~sunpy.net.solarnet.SOLARNETClient`
     - Ground-based high-resolution European solar observations (GREGOR, SST, THEMIS, VTT) via the SOLARNET VO.
     - `a.Time`, `a.Instrument`, `a.Wavelength`
   * - **XRSClient** / **GOESClient**
     - `~sunpy.net.dataretriever.XRSClient`
     - GOES X-ray Sensor (XRS) solar soft X-ray flux time series (Level 1, Level 2, netCDF & FITS archives).
     - `a.Time`, `a.Instrument('xrs')`
   * - **SUVIClient**
     - `~sunpy.net.dataretriever.SUVIClient`
     - GOES-R series (GOES-16, 17, 18, 19) Solar Ultraviolet Imager full-disk EUV imagery.
     - `a.Time`, `a.Instrument('suvi')`, `a.Wavelength`
   * - **EVEClient**
     - `~sunpy.net.dataretriever.EVEClient`
     - SDO/EVE Level 0CS Extreme Ultraviolet Variability Experiment irradiance data.
     - `a.Time`, `a.Instrument('eve')`, `a.Level`
   * - **GBMClient**
     - `~sunpy.net.dataretriever.GBMClient`
     - Fermi Gamma-Ray Burst Monitor (GBM) solar flare observations and lightcurves.
     - `a.Time`, `a.Instrument('gbm')`
   * - **LYRAClient**
     - `~sunpy.net.dataretriever.LYRAClient`
     - PROBA2/LYRA solar ultraviolet radiometer irradiance time series.
     - `a.Time`, `a.Instrument('lyra')`, `a.Level`
   * - **NoRHClient**
     - `~sunpy.net.dataretriever.NoRHClient`
     - Nobeyama RadioHeliograph averaged correlation radio flux time series (17 GHz & 34 GHz).
     - `a.Time`, `a.Instrument('norh')`, `a.Wavelength`
   * - **RHESSIClient**
     - `~sunpy.net.dataretriever.RHESSIClient`
     - RHESSI observing summary time series and flare catalogs.
     - `a.Time`, `a.Instrument('rhessi')`
   * - **GONGClient**
     - `~sunpy.net.dataretriever.GONGClient`
     - National Solar Observatory GONG synoptic magnetogram products.
     - `a.Time`, `a.Instrument('gong')`, `a.Physobs`
   * - **ADAPTClient**
     - `~sunpy.net.dataretriever.ADAPTClient`
     - NSO Air Force Data Assimilation Photospheric Flux Transport (ADAPT) magnetogram models.
     - `a.Time`, `a.Instrument('adapt')`
   * - **AIASynopsisClient**
     - `~sunpy.net.dataretriever.AIASynopsisClient`
     - AIA synoptic full-disk images from JSOC (accessible without email registration).
     - `a.Time`, `a.Instrument('aia')`, `a.Wavelength`
   * - **NOAAIndicesClient**
     - `~sunpy.net.dataretriever.NOAAIndicesClient`
     - NOAA Space Weather Prediction Center (SWPC) solar cycle indices (sunspot numbers, 10.7 cm radio flux).
     - `a.Time`, `a.Instrument('noaa-indices')`
   * - **NOAAPredictClient**
     - `~sunpy.net.dataretriever.NOAAPredictClient`
     - NOAA SWPC predicted solar cycle sunspot numbers and radio flux.
     - `a.Time`, `a.Instrument('noaa-predict')`
   * - **SRSClient**
     - `~sunpy.net.dataretriever.SRSClient`
     - NOAA SWPC Solar Region Summary (SRS) reports describing active regions and sunspot groups.
     - `a.Time`, `a.Instrument('srs-table')`
   * - **HEKClient**
     - `~sunpy.net.hek.HEKClient`
     - Heliophysics Event Knowledgebase (HEK) feature and event reports (flares, coronal holes, CMEs).
     - `a.Time`, `a.hek.EventType`
   * - **HECClient**
     - `~sunpy.net.helio.HECClient`
     - HELIO Web Services Event Catalogue search.
     - `a.Time`, `a.helio.TableName`

Detailed Client Guides
======================

Virtual Solar Observatory (VSO)
-------------------------------
The `~sunpy.net.vso.VSOClient` provides unified querying and downloading across dozens of international solar physics repositories.
It is often the primary source for full-disk and regional solar images.

* **Key Missions & Observatories Supported:**
  * SDO (AIA, HMI, EVE)
  * SOHO (EIT, LASCO, MDI, CDS, SUMER)
  * STEREO (EUVI, COR1, COR2, HI1, HI2)
  * Hinode (EIS, XRT, SOT)
  * Ground-based networks (BBSO, Kanzelhöhe, MLSO, TON)
* **Example Query:**

  .. code-block:: python

      >>> from sunpy.net import Fido, attrs as a
      >>> import astropy.units as u

      >>> results = Fido.search(  # doctest: +REMOTE_DATA
      ...     a.Time("2020-01-01 00:00", "2020-01-01 00:05"),
      ...     a.Instrument.aia,
      ...     a.Wavelength(171 * u.angstrom)
      ... )

Joint Science Operations Center (JSOC)
--------------------------------------
The `~sunpy.net.jsoc.JSOCClient` connects to the Stanford JSOC database hosting high-rate SDO data (AIA and HMI).
Unlike the VSO (which provides archived Level 1 files), JSOC enables querying data series, requesting custom cutouts, and exporting large quantities of data.

.. note::
   A registered email address with JSOC is required to export data. You can register your email address at `<http://jsoc.stanford.edu/ajax/register_email.html>`__.

* **Example Query:**

  .. code-block:: python

      >>> from sunpy.net import Fido, attrs as a

      >>> results = Fido.search(  # doctest: +REMOTE_DATA
      ...     a.Time("2014-01-01T00:00:00", "2014-01-01T00:01:00"),
      ...     a.jsoc.Series("hmi.m_45s"),
      ...     a.jsoc.Notify("user@example.com")
      ... )

For a complete guide on JSOC queries and staging requests, see :ref:`sunpy-tutorial-acquiring-data-jsoc`.

Solar Orbiter Archive (SOAR)
----------------------------
The `~sunpy.net.soar.SOARClient` interfaces with the European Space Agency's Solar Orbiter Archive.
It covers all science data levels (Level 0 through Level 3) across the entire mission payload.

* **Supported Payloads:**
  * Remote-sensing: EUI, PHI, SPICE, STIX, Metis, SoloHI
  * In-situ: MAG, EPD, RPW, SWA
* **Example Query:**

  .. code-block:: python

      >>> from sunpy.net import Fido, attrs as a

      >>> results = Fido.search(  # doctest: +REMOTE_DATA
      ...     a.Time("2021-02-01", "2021-02-02"),
      ...     a.Instrument("EUI"),
      ...     a.Level(1)
      ... )

Coordinated Data Analysis Web (CDAWeb)
--------------------------------------
The `~sunpy.net.cdaweb.CDAWEBClient` queries NASA's CDAWeb service, offering access to multi-mission heliophysics in-situ time series data, solar wind plasma parameters, and interplanetary magnetic field measurements.

* **Supported Missions:** Parker Solar Probe, Wind, ACE, STEREO, Voyager, MMS, and more.
* **Example Query:**

  .. code-block:: python

      >>> from sunpy.net import Fido, attrs as a

      >>> results = Fido.search(  # doctest: +REMOTE_DATA
      ...     a.Time("2020-06-01", "2020-06-02"),
      ...     a.cdaweb.Dataset("PSP_FLD_L2_MAG_RTN")
      ... )

Dedicated Mission Retrievers (``sunpy.net.dataretriever``)
----------------------------------------------------------
For many popular spacecraft and standard data feeds, ``sunpy`` provides lightweight, direct download clients in `sunpy.net.dataretriever`.
These clients fetch data directly from institutional FTP/HTTP archives without passing through the VSO.

* **GOES XRS (`~sunpy.net.dataretriever.XRSClient`):**
  Full mission coverage of GOES soft X-ray flux, used universally for solar flare classification.
* **GOES SUVI (`~sunpy.net.dataretriever.SUVIClient`):**
  Continuous full-disk extreme ultraviolet solar observations from the GOES-R series satellites.
* **PROBA-2 LYRA (`~sunpy.net.dataretriever.LYRAClient`):**
  High-cadence solar irradiance observations in four ultraviolet passbands.
* **Fermi GBM (`~sunpy.net.dataretriever.GBMClient`):**
  High-energy solar flare X-ray and gamma-ray lightcurves.
* **NOAA SWPC Indices (`~sunpy.net.dataretriever.NOAAIndicesClient`, `~sunpy.net.dataretriever.NOAAPredictClient`, `~sunpy.net.dataretriever.SRSClient`):**
  Daily and monthly solar cycle progression, sunspot counts, 10.7 cm radio fluxes, and active region summaries.

Heliophysics Event Knowledgebase (HEK)
--------------------------------------
The `~sunpy.net.hek.HEKClient` queries event catalogs rather than raw image files.
It allows users to find when and where specific solar events occurred (such as flares, coronal mass ejections, or filament eruptions).

* **Example Query:**

  .. code-block:: python

      >>> from sunpy.net import Fido, attrs as a

      >>> results = Fido.search(  # doctest: +REMOTE_DATA
      ...     a.Time("2011-06-07 06:00", "2011-06-07 08:00"),
      ...     a.hek.EventType("FL")
      ... )

For a detailed tutorial on filtering event results, see :ref:`sunpy-tutorial-acquiring-data-hek`.

Tips for Querying Clients with Fido
===================================

1. **Check Supported Values for Attributes:**
   You can print any attribute class to see all known valid parameters recognized by ``sunpy``:

   .. code-block:: python

       >>> print(a.Instrument)
       >>> print(a.Source)

2. **Targeting a Specific Client:**
   If a general query like ``Fido.search(a.Time(...), a.Instrument.aia)`` returns results from both VSO and another provider and you wish to restrict results to a specific client, you can filter by provider or use client-specific attributes.

3. **Downloading Your Search Results:**
   Once you have inspected your search results, download them using `~sunpy.net.Fido.fetch`:

   .. code-block:: python

       >>> files = Fido.fetch(results)  # doctest: +SKIP
