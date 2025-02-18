.. _sunpy-how-to-plot-atmos-model:

*****************************
Plot a solar atmosphere model
*****************************

Here we will demonstrate how to load and visualize one of the models in `sunpy.sun.models`, in this case :cite:t:`avrett_loeser_2008` solar chromosphere model.

.. plot::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import sunpy.sun.models as sun_models

    >>> data = models.avrett_loeser_2008

    >>> x_param = "h"
    >>> y_param = ["T","n_e"]

    >>> x_values = data[x_param]
    >>> y_units = {param: f" ({data[param].unit})" if hasattr(data[param], "unit") else "" for param in y_params}

    >>> x_unit = f" ({x_values.unit})" if hasattr(x_values, "unit") else ""
    >>>y_units = {param: f" ({data[param].unit})" if hasattr(data[param], "unit") else "" for param in y_params}

    >>>plt.figure(figsize=(8, 5))

    >>>plt.plot(x_values, y_values["T"], marker="o", linestyle="-", label=f"Temperature {y_units['T']}")

    >>>plt.plot(x_values, y_values["n_e"], marker="s", linestyle="--", label=f"Density {y_units['n_e']}")

    >>>plt.yscale("log")
    >>>plt.xscale("log")


    >>>plt.xlabel(f"{x_param}{x_unit}")
    >>>plt.ylabel("Temperature & Density (log scale)")
    >>>plt.title("Solar Chromosphere Model: Temperature & Density vs Height")

    >>>plt.legend()
    >>>plt.grid(True, which="both", linestyle="--", linewidth=0.5)

    >>>plt.show()
