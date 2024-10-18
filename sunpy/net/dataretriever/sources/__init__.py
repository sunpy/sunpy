"""
Sources for the SunPy DataRetriever.

This module contains client classes for various data sources, allowing users
to easily query and retrieve solar data from multiple databases and services.
"""

from .adapt import *
from .aia_synoptic import *
from .eve import *
from .fermi_gbm import *
from .goes import *
from .gong import *
from .lyra import *
from .noaa import *
from .norh import *
from .rhessi import *

# List of all client classes available in this module
__all__ = [
    "ADAPTClient",
    "AIASynopticClient",
    "EVEClient",
    "FermiGBMClient",
    "GOESClient",
    "GONGClient",
    "LYRAClient",
    "NOAAIndicesClient",
    "NOAAPredictClient",
    "NOAASRSClient",
    "NoRHClient",
    "RHESSIClient",
]
