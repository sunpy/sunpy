"""
===============
sunpy Constants
===============

The example shows all of the solar physics specific constants provided by sunpy.
"""

from sunpy.sun import constants as con

##############################################################################
# All constants are stored in a dictionary. A list of all available keys
# which describe each constant can be had with the following command.

print(con.constants.keys())

##############################################################################
# The following command will display all 34 constants.

print(con.print_all())
