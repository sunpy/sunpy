# coding: utf-8
"""
===============
SunPy Constants
===============

The example shows all of the constants provided by sunpy.
"""

from sunpy.sun import constants as con

##############################################################################
# All constants are stored in a dictionary. A list of all available keys
# which describe the constant can be had with the following command.
con.constants.keys()

##############################################################################
# The following command will display all constants as well as their values
# in an astropy `Table <http://docs.astropy.org/en/stable/table/index.html>`_
print(con.print_all())
