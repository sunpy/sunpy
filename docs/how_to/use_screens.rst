.. _sunpy-how-to-use_screens:

*********************************
Specify any use of screens (Plot Windows)
*********************************

When working with SunPy, you often visualize solar data using maps, time series, or spectra. These are displayed as **plot windows**, or "screens", generated using Matplotlib.

These screens are *temporary visualizations* created at runtime. They are helpful for:
- Viewing solar features (e.g., sunspots, flares, active regions)
- Inspecting data from different instruments or wavelengths
- Debugging analysis steps visually

Basic Example
-------------

Here's how to load and display a solar image using SunPy:

.. code-block:: python

    import sunpy.map
    import sunpy.data.sample
    import matplotlib.pyplot as plt

    # Load a sample AIA image
    smap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

    # Display the image
    smap.plot()
    plt.colorbar()
    plt.title("Solar Image in AIA 171 Å")
    plt.show()

The above code will open a plot window if you're using a desktop Python environment (e.g., VS Code, Anaconda, or running from a script). If you're in a Jupyter notebook, the plot will appear inline.

Saving the Plot (Optional)
--------------------------

Screens are not saved automatically. If you'd like to keep the plot for reports or later use, you can save it as an image file:

.. code-block:: python

    plt.savefig("aia171_plot.png")

.. note::
    Screens (plot windows) will disappear once closed or when the session ends unless saved using ``plt.savefig()``.

Common Pitfalls
---------------

- **Closing the screen window before saving** will result in lost work.
- **Interactive tools like Jupyter** do not persist plots between sessions unless saved.
- Always call ``plt.show()`` to display the plot when running scripts outside of notebooks.

