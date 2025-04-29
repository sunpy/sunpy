"""
=====================================
Plotting the Solar Atmosphere Model
=====================================

This example demonstrates how to load and visualize the
Avrett & Loeser (2008) solar atmosphere model using SunPy.
"""

import matplotlib.pyplot as plt

import sunpy.sun.models as sun_models

###############################################################################
# Load the atmosphere model from SunPy's `sun.models` module.

data = sun_models.chromosphere_avrett_loeser_2008

###############################################################################
# Define the X and Y parameters to plot.
# We'll use height (`h`) on the X-axis and temperature (`T`) and
# electron density (`n_e`) on the Y-axes.

x_param = "h"
y_param = ["T", "n_e"]

x_unit = f" ({data[x_param].unit})" if hasattr(data[x_param], "unit") else ""
y_units = {
    param: f" ({data[param].unit})" if hasattr(data[param], "unit") else ""
    for param in y_param
}

x_values = data[x_param]
y_values = {param: data[param] for param in y_param}

###############################################################################
# Create a dual-axis plot for temperature and electron density.
# Temperature will be shown on the left Y-axis and electron density on the right.

fig, ax1 = plt.subplots(figsize=(8, 5))

ax1.plot(x_values, y_values["T"], color="tab:red", marker="o", linestyle="-")
ax1.set_xlabel(f"{x_param}{x_unit}")
ax1.set_ylabel(f"Temperature{y_units['T']}", color="tab:red")
ax1.tick_params(axis="y", labelcolor="tab:red")

ax2 = ax1.twinx()
ax2.plot(x_values, y_values["n_e"], color="tab:blue", marker="s", linestyle="--")
ax2.set_ylabel(f"Electron Density{y_units['n_e']}", color="tab:blue")
ax2.tick_params(axis="y", labelcolor="tab:blue")

###############################################################################
# Apply logarithmic scaling to both axes for better visualization of the data.

ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_yscale("log")

###############################################################################
# Add a title and grid for clarity.

plt.title("Avrett & Loeser (2008) Solar Atmosphere Model")
ax1.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.show()
