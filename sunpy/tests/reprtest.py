from sunpy.net.vso.attrs import Wave
import pytest

def test_repr():
    """Calls the __repr__ method of class vso.attrs.Wave"""
    wav = Wave(12, 16)
    print(wav)

test_repr()
