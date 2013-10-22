from __future__ import absolute_import

__all__ = ['test']

import sunpy
import pytest

def test():
    pytest.main(sunpy.__path__)

