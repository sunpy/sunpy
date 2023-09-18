``sunpy.map.extract_along_coord()`` has been removed.
Instead, use :func:`~sunpy.map.pixelate_coord_path`, and then pass its output to :func:`~sunpy.map.sample_at_coords`.
``pixelate_coord_path`` uses a different line algorithm by default, but you can specify ``bresenham=True`` as an argument to use the same line algorithm as ``extract_along_coord``.
