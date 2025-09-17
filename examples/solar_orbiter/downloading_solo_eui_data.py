"""
========================================
Downloading Solar Orbiter EUI Data
========================================

This example shows how to download Solar Orbiter Extreme Ultraviolet Imager (EUI)
data using the `sunpy-soar <https://pypi.org/project/sunpy-soar/>`__ package.

The sunpy-soar package enables searching and downloading data from the
Solar Orbiter Archive (SOAR) through the `~sunpy.net.Fido` interface.
"""
import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, we need to define our search parameters. Let's search for EUI images
# taken around the time of a solar event. We'll look for Full Sun Imager (FSI)
# data in the 174 Angstrom channel.

time_range = a.Time('2022-03-30T10:00:00', '2022-03-30T11:00:00')
instrument = a.Instrument('EUI')
level = a.Level(2)  # Level 2 processed data
product_type = a.soar.ProductType('EUI-FSI174-IMAGE')

###############################################################################
# Now we can search for the data using Fido. The sunpy-soar package automatically
# registers with Fido when imported, so we can search the Solar Orbiter Archive
# directly.

try:
    import sunpy_soar  # noqa: F401
    result = Fido.search(time_range, instrument, level, product_type)
    print(result)

    ###############################################################################
    # Let's download the first file from our search results.

    if len(result) > 0:
        files = Fido.fetch(result[0])
        print(f"Downloaded: {files}")

        ###############################################################################
        # Now let's load the downloaded data and create a Map.

        eui_map = sunpy.map.Map(files[0])
        print(f"Created map: {eui_map}")
        print(f"Observation time: {eui_map.date}")
        print(f"Instrument: {eui_map.instrument}")
        print(f"Observatory: {eui_map.observatory}")
        print(f"Wavelength: {eui_map.wavelength}")

    else:
        print("No data found for the specified search criteria.")

except ImportError:
    print("This example requires the sunpy-soar package to be installed.")
    print("You can install it with: pip install sunpy-soar")

    ###############################################################################
    # Alternative: Show how to install and import sunpy-soar

    print("\nAlternatively, you can install sunpy-soar using conda:")
    print("conda install -c conda-forge sunpy-soar")

    ###############################################################################
    # For demonstration purposes, let's show what the import structure would look like:

    print("\nOnce installed, you would import sunpy-soar like this:")
    print("import sunpy_soar")
    print("from sunpy.net import attrs as a")
    print("result = Fido.search(a.Time(...), a.Instrument('EUI'), a.soar.ProductType(...))")
