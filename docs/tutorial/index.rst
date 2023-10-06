.. _sunpy-tutorial-index:

******************
The sunpy tutorial
******************

Welcome to the introductory tutorial for the ``sunpy`` core package.
``sunpy`` is a community-developed, free and open-source solar data analysis environment.
It is meant to provide the core functionality and tools to analyze solar data with Python.

**Who this is for**

This tutorial assumes you know how to run code in Python, but doesn't make any assumptions beyond that.
``sunpy`` depends heavily on other packages in the scientific Python ecosystem, including Astropy, NumPy, and Matplotlib.
This tutorial explains the basics of these packages needed to use ``sunpy``, but is not a comprehensive tutorial for these packages and references their own introductory tutorials where relevant.

**How to use this tutorial**

This tutorial is designed to be read from start to finish.
You can either read it without running any code, or copy and paste the code snippets as you read.
Each chapter of the tutorial provides a self-contained set of codes.

Later chapters build on the concepts introduced in earlier chapters.
The one exception is the two chapters that explore Maps and TimeSeries - these are independent, but depend on all earlier chapters.

.. toctree::
   :maxdepth: 2

   installation
   units
   time
   coordinates
   acquiring_data/index
   maps
   timeseries
