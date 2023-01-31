#!/bin/env python

import asyncio

import sunpy
from sunpy.net.cdaweb.helpers import _update_cdaweb_dataset_data

print(f"Updating the attrs json files using sunpy {sunpy.__version__}...")

print("Updating CDAWeb json...\nThis may take some time...")
asyncio.run(_update_cdaweb_dataset_data())

print("Done. Don't forget to update the doctests.")
