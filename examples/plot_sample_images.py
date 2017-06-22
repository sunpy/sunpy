"""
=============
Sample Images
=============

This example shows off all of the sample images provided by SunPy.
"""
import matplotlib.pyplot as plt

from sunpy.data.sample import file_dict
from sunpy.map import Map

###############################################################################
# We create a loop through all of the sample shortcuts and picks out the images
# that can be opened by Map. It then plots each one.
for key in file_dict:
    if key.count('IMAGE'):
        try:
            m = Map(file_dict[key])
            plt.figure()
            m.plot()
            plt.savefig('{0}.png'.format(key))
            #plt.show()
        except:
            pass
