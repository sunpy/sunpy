import io
import json
import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

import astropy

__all__ = ["hash_figure", "verify_figure_hash"]

ft2_version = f"{mpl.ft2font.__freetype_version__.replace('.', '')}"
mpl_version = "dev" if "+" in mpl.__version__ else mpl.__version__.replace('.', '')
astropy_version = "dev" if "dev" in astropy.__version__ else astropy.__version__.replace('.', '')
HASH_LIBRARY_NAME = f"figure_hashes_mpl_{mpl_version}_ft_{ft2_version}_astropy_{astropy_version}.json"
HASH_LIBRARY_FILE = Path(__file__).parent / HASH_LIBRARY_NAME

# Load the hash library if it exists
try:
    with open(HASH_LIBRARY_FILE) as infile:
        print(f"Figure tests are using hash library: {HASH_LIBRARY_NAME}")
        hash_library = json.load(infile)
except OSError:
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
