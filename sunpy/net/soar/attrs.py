"""
This file defines the net attributes that can be used to search the SOAR.
"""
from ._attrs import SOOP, Distance, Product

__all__ = ["SOOP", "Distance", "Product"]

# Trick the docs into thinking these attrs are defined in here.
for _a in (SOOP, Distance, Product):
    _a.__module__ = __name__
