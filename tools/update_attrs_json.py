#!/bin/env python

import sunpy
from sunpy.net.jsoc import JSOCClient
from sunpy.net.vso import VSOClient

print(f"Updating the attrs json files using sunpy {sunpy.__version__}...")


print("Updating VSO json...")

VSOClient.create_parse_vso_values()


print("Updating JSOC json...\nThis may take some time...")

JSOCClient.create_parse_jsoc_values()

print("Done. Don't forget to update the doctests.")
