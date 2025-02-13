.. _sunpy-plot-chromosphere-model:

***********************************
Plotting the Chromosphere Model
***********************************

This example demonstrates how to load and visualize the
Avrett & Loeser (2008) solar chromosphere model data using Astropy.

You can specify custom x and y parameters and optionally apply
a log scale to the y-axis.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import sunpy.sun.models as models

    >>> data = models.chromosphere_avrett_loeser_2008

    >>> x_param = "h"
    >>> y_param = "T"
    >>> log = True
    >>> color = "b"

    >>> if x_param not in data.colnames or y_param not in data.colnames:
    >>>     raise ValueError(f"Invalid parameters. Available options: {list(data.colnames)}")

    >>> x_values = data[x_param]
    >>> y_values = data[y_param]

    >>> y_label = y_param
    >>> if log:
    >>>     y_values = np.log10(y_values.value)
    >>>     y_label = f"log10({y_param})"

    >>> x_unit = f"({x_values.unit})" if hasattr(x_values, "unit") else ""
    >>> y_unit = f"({data[y_param].unit})" if hasattr(data[y_param], "unit") else ""

    >>> plt.figure(figsize=(8, 5))
    >>> plt.plot(x_values, y_values, marker="o", linestyle="-", color=color, label=f"{y_label} vs {x_param}")
    >>> plt.xlabel(f"{x_param} {x_unit}")
    >>> plt.ylabel(f"{y_label} {y_unit}")
    >>> plt.title(f"Solar Chromosphere Model: {y_label} vs {x_param}")
    >>> plt.legend()
    >>> plt.grid(True)
    >>> plt.show()
