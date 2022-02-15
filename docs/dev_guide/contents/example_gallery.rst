.. _example_gallery:

***************
Example Gallery
***************

The purpose of the page is to describe the contribution guidelines for the `sunpy Example Gallery <https://docs.sunpy.org/en/stable/generated/gallery/index.html>`_.

All potential contributors to the sunpy Example Gallery should read and abide by the following guidelines.

.. note:: We have an example template located at ``examples/example_template/example_template.py``.

Contribution Guidelines
=======================

* The title of the example should be short yet descriptive and emphasize the goal of the example.
  Try to make the title appeal to a broad audience and avoid referencing a specific instrument, catalog, or anything wavelength dependent.

* Each example should begin with a paragraph that gives a brief overview of the entire example, including relevant astronomy concepts, and motivates the described functionality.

* The examples must be compatible with the versions supported by the last major release of the sunpy core package (i.e., Python >= 3.7).

* All the examples must be fully PEP8 compliant, we recommend using one of the many PEP8 linters that are available (autopep8, flake8 as some examples).

* Wherever possible, the examples should include linked references with links pointing to the appropriate `DOI <https://zenodo.org/record/2551710>`_ or `ADS <https://ui.adsabs.harvard.edu/>`_ entry.

* The example should include links to relevant documentation pages.

* Each example should, where possible, include at least one image, map, or plot to use as the icon in the example gallery.

* The examples should avoid using acronyms without defining them first (e.g. Virtual Solar Observatory, or VSO).
  Similarly complex jargon should be avoided unless clearly explained.

* There should be a good variety of examples for each section (simple and more complex to cater for different levels).

* When creating a plot, particularly of a map, the example should follow these rules of thumb to minimize verbosity and maintain consistency with the other gallery examples:

  * Do not use `~sunpy.map.GenericMap.peek` in examples. Instead, use `~sunpy.map.GenericMap.plot`.

  * Always create a figure instance using ``plt.figure()`` prior to creating a plot.

  * Only create an axes instance if it is explicitly needed later in the example
    (e.g. when overplotting a coordinate on a map using ``ax.plot_coord``).

  * If an axes instance is created, it should be explicitly passed as a keyword argument wherever possible
    (e.g. in `~sunpy.map.GenericMap.plot` or `~sunpy.map.GenericMap.draw_grid`).

  * After each figure is created, call ``plt.show()``. While not explicitly needed for the gallery to render,
    this ensures that, when run as a script, the gallery will display each figure sequentially.
