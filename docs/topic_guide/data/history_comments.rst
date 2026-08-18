.. _sunpy-topic-guide-history-comments:

**********************************************************
How "HISTORY" and "COMMENT" FITS Keys are Handled in sunpy
**********************************************************

In FITS files, there are often multiple entries for both the "HISTORY" and "COMMENT" keys.
For example, when applying a prep routine to an image, a "HISTORY" entry may be added to the FITS header for every step in the prep pipeline.
Because the metadata associated with each `~sunpy.map.GenericMap` acts like a dictionary, where each key must be unique, these repeated "HISTORY" and "COMMENT" keys cannot be represented as separate entries in `~sunpy.util.MetaDict`.
Thus, when a FITS file with multiple "HISTORY" keys is read into a Map object, all the values corresponding to "HISTORY" are joined together with newline characters ``\n`` and stored as a single "history" entry in `~sunpy.util.MetaDict`.
When writing the resulting Map to a FITS file, "history" is split along the ``\n`` characters and each entry is written to a separate "HISTORY" key in the resulting FITS header.
The same is true for "COMMENT" keys.

`See this pull request for additional details. <https://github.com/sunpy/sunpy/pull/6911>`__
