"""
This file defines the net attributes that can be used to search the SOAR.
"""

from ._attrs import SOOP, Category, Distance, Product, Sensor

__all__ = ["SOOP", "Distance", "Product", "Sensor", "Category"]

# Trick the docs into thinking these attrs are defined in here.
for _a in (SOOP, Distance, Product, Sensor, Category):
    _a.__module__ = __name__
