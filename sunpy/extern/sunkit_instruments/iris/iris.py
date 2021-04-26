"""
This package provides Interface Region Imaging Spectrometer (IRIS) instrument
routines.

.. note::

    More comprehensive IRIS tools are now being developed by the
    `IRIS instrument team. <https://gitlab.com/LMSAL_HUB/iris_hub>`__
"""
import sunpy.io
import sunpy.map
import sunpy.time

__all__ = ["SJI_to_sequence"]


def SJI_to_sequence(filename, start=0, stop=None, hdu=0):
    """
    Read a SJI file and return a `sunpy.map.MapSequence`.

    .. warning::
        This function is a very early beta and is not stable. Further work is
        on going to improve SunPy IRIS support.

    Parameters
    ----------
    filename: `str`
        File to read.
    start: `int`, optional
        Temporal axis index to create `~sunpy.map.MapSequence` from.
        Defaults to 0, which will start from the begining.
    stop: `int`, optional
        Temporal index to stop `~sunpy.map.MapSequence` at.
        Defaults to `None`, which will use the entire index.
    hdu: `int`, optional
        The hdu index to use, defaults to 0.

    Returns
    -------
    `~sunpy.map.MapSequence`
        A map sequence of the SJI data.
    """

    hdus = sunpy.io.read_file(filename)
    # Get the time delta
    time_range = sunpy.time.TimeRange(hdus[hdu][1]["STARTOBS"], hdus[hdu][1]["ENDOBS"])
    splits = time_range.split(hdus[hdu][0].shape[0])

    if not stop:
        stop = len(splits)

    headers = [hdus[hdu][1]] * (stop - start)
    datas = hdus[hdu][0][start:stop]

    # Make the cube:
    iris_cube = sunpy.map.Map(list(zip(datas, headers)), sequence=True)
    # Set the date/time

    for i, m in enumerate(iris_cube):
        m.meta["DATE-OBS"] = splits[i].center.isot

    return iris_cube
