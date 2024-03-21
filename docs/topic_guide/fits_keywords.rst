.. _sunpy-topic-guide-fits-keywords:

************************************
Creating and validating FITS Headers
************************************

Sometimes it is necessary to inspect, validate, or create metadata for a FITS file.
There are a number of standards and requirements associated with the keywords used in FITS files.
When opening reading a FITS file, sunpy makes use of these standards to properly interpret the FITS metadata.

The SOLARNET Metadata Recommendations for Solar Observations (Version 2.0, 14. August 2023) `https://arxiv.org/pdf/2011.12139.pdf <https://arxiv.org/pdf/2011.12139.pdf>`_ provides descriptions and recommendations for many keywords.
A listing of FITS keywords and their descriptions are provided below.
The `required` column describes whether the keyword is required by SOLARNET.
The `data_type` column provides a description of the expected type of the keyword value.
Note that deprecated keywords (such as DATE-OBS or EXPTIME) are explicitely not included in this list.


.. csv-table:: Table 1-1: FITS keywords
   :file: ../generated/fits_schema.csv
   :widths: 30, 70, 30, 30, 30
   :header-rows: 1


In order to support creating FITS headers with the appropriate keywords, the
`sunpy.io.meta.fits_meta` module provides the `~sunpy.io.meta.fits_meta.SolarnetHeader` class.
It is a subclass of `~astropy.io.fits.Header`.
It provides additional functionality like the ability to validate against the solarnet schema as well as support an additional custom schema.

To create a blank header with only the required keywords

.. code-block:: python

    from sunpy.io.meta import fits_meta
    header = fits_meta.SolarnetHeader()
    header['OBSRVTRY'] = 'myawesomemission'

You can then inspect the keywords and fill in the values as needed.
After you've filled all of your values, you can validate the header with

.. code-block:: python

    header.validate()

This will provide warnings for each issue that it finds.
It checks for a number of potential issues such as whether the required keywords 
have non-blank values, that deprecrated keywords are not used, and that the values 
can be interpreted as the correct types (e.g. dates, quantities, int, floats).
For a complete listing of the issues checked see `~sunpy.io.meta.fits_meta.SolarnetHeader.validate()`.

You can define your own schema by creating a yaml file with the following structure for each keywords.

.. code-block:: yaml

    SOLARNET:
        description: Fully SOLARNET-compliant=1.0, partially=0.5
        human_readable: Solarnet compatibility
        required: true
        data_type: float

You can add an additional optional field `allowed_values: [1.0, 0.5, 0]` if you want to limit the possible options.