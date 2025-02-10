"""
=========================================
Plotting the Chromosphere Model
=========================================

This example demonstrates how to load and visualize the
Avrett & Loeser (2008) solar chromosphere model data using Astropy.

You can specify custom x and y parameters, and optionally apply
a log scale to the y-axis.
"""

import matplotlib.pyplot as plt
import numpy as np

import sunpy.sun.models as models

# Load the chromosphere data
data = models.get_chromosphere_data

# Specify x and y parameters
x_param = "h"   # Height
y_param = "T"   # Temperature
log = True      # Apply log scale to y-axis

# Validate parameters
if x_param not in data.colnames or y_param not in data.colnames:
    raise ValueError(f"Invalid parameters. Available options: {list(data.colnames)}")

# Extract values
x_values = data[x_param]
y_values = data[y_param]

# Apply log scale if needed
if log:
    y_values = np.log10(y_values)
    y_label = f"log10({y_param})"
else:
    y_label = y_param

# Plot
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, marker="o", linestyle="-", color="b", label=f"{y_label} vs {x_param}")

# Handle units
x_unit = f"({data[x_param].unit})" if data[x_param].unit else ""
y_unit = f"({data[y_param].unit})" if data[y_param].unit else ""

plt.xlabel(f"{x_param} {x_unit}")
plt.ylabel(f"{y_label} {y_unit}")
plt.title(f"Solar Chromosphere Model: {y_label} vs {x_param}")
plt.legend()
plt.grid(True)
plt.show()

# Example: Custom plot with a linear scale
x_param = "m"   # Column name for x-axis
y_param = "n_e" # Column name for y-axis
log = False     # No log scale

# Extract new values
x_values = data[x_param]
y_values = data[y_param]

# Update y-axis label
y_label = y_param

# Plot again with new parameters
plt.figure(figsize=(8, 5))
plt.plot(x_values, y_values, marker="o", linestyle="-", color="r", label=f"{y_label} vs {x_param}")

# Handle units
x_unit = f"({data[x_param].unit})" if data[x_param].unit else ""
y_unit = f"({data[y_param].unit})" if data[y_param].unit else ""

plt.xlabel(f"{x_param} {x_unit}")
plt.ylabel(f"{y_label} {y_unit}")
plt.title(f"Solar Chromosphere Model: {y_label} vs {x_param}")
plt.legend()
plt.grid(True)
plt.show()
