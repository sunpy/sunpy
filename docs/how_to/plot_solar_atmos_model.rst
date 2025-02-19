.. _sunpy-how-to-plot-atmos-model:

*****************************
Plot a solar atmosphere model
*****************************

Here we will demonstrate how to load and visualize one of the models in `sunpy.sun.models`, in this case :cite:t:`avrett_loeser_2008` solar chromosphere model.

.. plot::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import sunpy.sun.models as sun_models

    >>> data = sun_models.avrett_loeser_2008

    >>> x_param = "h"
    >>> y_param = ["T","n_e"]

    >>> x_values = data[x_param]
    >>> y_units = {param: f" ({data[param].unit})" if hasattr(data[param], "unit") else "" for param in y_params}

    >>> fig, ax1 = plt.subplots(figsize=(8, 5))

    >>> ax1.plot(x_values, y_values["T"], marker="o", linestyle="-", color="tab:red", label=f"Temperature {y_units['T']}")
    >>> ax1.set_xlabel(f"{x_param}{x_unit}")
    >>> ax1.set_ylabel(f"Temperature {y_units['T']}", color="tab:red")
    >>> ax1.tick_params(axis="y", labelcolor="tab:red")

    >>> ax2 = ax1.twinx()
    >>> ax2.plot(x_values, y_values["n_e"], marker="s", linestyle="--", color="tab:blue", label=f"Density {y_units['n_e']}")
    >>> ax2.set_ylabel(f"Electron Density {y_units['n_e']}", color="tab:blue")
    >>> ax2.tick_params(axis="y", labelcolor="tab:blue")

    >>> ax1.set_xscale("log")
    >>> ax1.set_yscale("log")
    >>> ax2.set_yscale("log")

    >>> plt.title("Avrett & Loeser (2008) Model")
    >>> ax1.grid(True, which="both", linestyle="--", linewidth=0.5)

    >>> plt.show()
