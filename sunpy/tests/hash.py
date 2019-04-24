import io
import os
import json
import hashlib
from sys import version_info

import matplotlib.pyplot as plt

__all__ = ["hash_figure", "verify_figure_hash"]

HASH_LIBRARY_NAME = 'figure_hashes_py{0}{1}.json'.format(version_info.major, version_info.minor)
HASH_LIBRARY_FILE = os.path.join(os.path.dirname(__file__), HASH_LIBRARY_NAME)

# Load the hash library if it exists
try:
    with open(HASH_LIBRARY_FILE) as infile:
        hash_library = json.load(infile)
except IOError:
    hash_library = {}


def hash_figure(figure=None, out_stream=None):
    """
    For a `matplotlib.figure.Figure`, returns the SHA256 hash as a hexadecimal
    string.

    Parameters
    ----------
    figure : `matplotlib.figure.Figure`, optional
        Defaults to `None` and the current figure is used.
    out_stream : I/O stream (e.g., an open file), optional
        Defaults to `None`. If not, write a PNG of the figure to the stream.

    Returns
    -------
    `str`
        The SHA256 hash in hexadecimal representation.
    """

    if figure is None:
        figure = plt.gcf()

    if out_stream is None:
        imgdata = io.BytesIO()
    else:
        imgdata = out_stream

    figure.savefig(imgdata, format='png')

    out = _hash_file(imgdata)
    if out_stream is None:
        imgdata.close()
    return out


def _hash_file(in_stream):
    """
    Hashes an already opened file.
    """
    in_stream.seek(0)
    buf = in_stream.read()
    hasher = hashlib.sha256()
    hasher.update(buf)
    return hasher.hexdigest()


def verify_figure_hash(name, figure=None):
    """
    Verifies whether a figure has the same hash as the named hash in the
    current hash library. If the hash library does not contain the specified
    name, the hash is added to the library.

    Parameters
    ----------
    name : `str`
        The identifier for the hash in the hash library.
    figure : `matplotlib.figure.Figure`, optional
        Defaults to `None`. if not, the current figure is used.

    Returns
    -------
    `bool`
        `False` if the figure's hash does not match the named hash, otherwise `True`.
    """
    if name not in hash_library:
        hash_library[name] = hash_figure(figure)
        return True
    return hash_library[name] == hash_figure(figure)
