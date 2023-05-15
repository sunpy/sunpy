.. _example_gallery:

***************
Example Gallery
***************

The purpose of the page is to describe the contribution guidelines for the `sunpy Example Gallery <https://docs.sunpy.org/en/stable/generated/gallery/index.html>`__.

All potential contributors to the ``sunpy`` example gallery should read and abide by the following guidelines.

Contribution Guidelines
=======================

* The title of the example should be short yet descriptive and emphasize the goal of the example.
  Try to make the title appeal to a broad audience and avoid referencing a specific instrument, catalog, or anything wavelength dependent.

* Each example should begin with a paragraph that gives a brief overview of the entire example, including relevant astronomy concepts, and motivates the described functionality.

* The examples must be compatible with the versions supported by the last major release of the ``sunpy`` core package.

* All the examples must be fully PEP8 compliant, the ``pre-commit`` hook should be used to ensure this.

* Wherever possible, the examples should include linked references with links pointing to the appropriate `DOI <https://zenodo.org/record/2551710>`__ or `ADS <https://ui.adsabs.harvard.edu/>`__ entry.

* The example should include links (URL or `sphinx intersphinx <https://coderefinery.github.io/sphinx-lesson/intersphinx/>`__) to relevant documentation pages.

* Each example should, where possible, include at least one image, map, or plot to use as the icon in the example gallery.

* The examples should avoid using acronyms without defining them first (e.g. Virtual Solar Observatory, or VSO).
  Similarly complex jargon should be avoided unless clearly explained.

* There should be a good variety of examples for each section (simple and more complex to cater for different levels).

* When creating a plot, particularly of a map, the example should follow these rules of thumb to minimize verbosity and maintain consistency with the other gallery examples:

  * Do not use `~sunpy.map.GenericMap.peek` in examples. Instead, use `~sunpy.map.GenericMap.plot`.

  * Always create a figure instance using ``plt.figure()`` prior to creating a plot.

  * Always create an Axes instance, and where possible use this to modify the plot instead of using the ``pyplot`` interface (for example use ``ax.set_xlabel()``` instead of ``plt.xlabel()``).

  * If an axes instance is created, it should be explicitly passed as a keyword argument wherever possible (e.g. in `~sunpy.map.GenericMap.plot` or `~sunpy.map.GenericMap.draw_grid`).

  * At the end of the example,  call ``plt.show()``.
    While not explicitly needed for the gallery to render, this ensures that when run as a script, the plots will appear.

  * If you need to use ``astropy.visualization.{quantity_support, time_support}``, import these functions at the top of the example, and call them directly before the first plot that needs them.

We recommend checking the other examples in the gallery to see how they are structured and what they contain.
