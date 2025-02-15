.. _sunpy-plot-chromosphere-model:

***********************************
Plotting the Chromosphere Model
***********************************

This example demonstrates how to load and visualize one of the models in `sunpy.sun.models`, in this case :cite:t:`avrett_loeser_2008` solar chromosphere model.
Which we will just plot a range of parameters (temperature and density as a function of height).

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import sunpy.sun.models as sun_models

    >>> data = models.avrett_loeser_2008

    >>> x_param = "h"
    >>> y_param = "T"

    >>> x_values = data[x_param]
    >>> y_values = data[y_param]

    >>> x_unit = f"({x_values.unit})" if hasattr(x_values, "unit") else ""
    >>> y_unit = f"({data[y_param].unit})" if hasattr(data[y_param], "unit") else ""

    >>> plt.figure(figsize=(8, 5))
    >>> plt.plot(x_values, y_values, marker="o", linestyle="-", label=f"{y_param} vs {x_param}")
    >>> plt.xlabel(f"{x_param} {x_unit}")
    >>> plt.ylabel(f"{y_param} {y_unit}")
    >>> plt.title(f"Solar Chromosphere Model: {y_param} vs {x_param}")
    >>> plt.legend()
    >>> plt.grid(True)
    >>> plt.show()
