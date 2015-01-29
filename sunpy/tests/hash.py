import os
import io
import hashlib
import json

import matplotlib.pyplot as plt
import pytest

HASH_LIBRARY_NAME = 'figure_hashes.json'

# Load the hash library if it exists
try:
    with open(os.path.join(os.path.dirname(__file__), HASH_LIBRARY_NAME)) as infile:
        hash_library = json.load(infile)
except IOError:
    hash_library = {}

def hash_figure(figure=None):
    """
    For a matplotlib.figure.Figure, returns the SHA256 hash as a hexadecimal string.

    Parameters
    ----------
    figure : matplotlib.figure.Figure
        If None is specified, the current figure is used (as determined by matplotlib.pyplot.gcf())

    Returns
    -------
    out : string
        The SHA256 hash in hexadecimal representation
    """

    if figure is None:
        figure = plt.gcf()

    imgdata = io.BytesIO()
    figure.savefig(imgdata, format='png')

    imgdata.seek(0)
    buf = imgdata.read()
    imgdata.close()

    hasher = hashlib.sha256()
    hasher.update(buf)
    return hasher.hexdigest()

def verify_figure_hash(name, figure=None):
    """
    Verifies whether a figure has the same hash as the named hash in the current hash library.
    If the hash library does not contain the specified name, the hash is added to the library.

    Parameters
    ----------
    name : string
        The identifier for the hash in the hash library
    figure : matplotlib.figure.Figure
        If None is specified, the current figure is used (as determined by matplotlib.pyplot.gcf())

    Returns
    -------
    out : bool
        False if the figure's hash does not match the named hash, otherwise True
    """
    if name not in hash_library:
        hash_library[name] = hash_figure(figure)
        return True
    return hash_library[name] == hash_figure(figure)

def assert_figure_hash(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the returned figure,
    with the name of the test function as the hash identifier in the library.
    """
    @pytest.mark.figure
    def wrapper(*args, **kwargs):
        name = test_function.func_name
        hash = hash_figure(test_function(*args, **kwargs))
        if name not in hash_library:
            hash_library[name] = hash
            pytest.fail("Hash not present: {0}".format(name))
        else:
            assert hash_library[name] == hash
    return wrapper
