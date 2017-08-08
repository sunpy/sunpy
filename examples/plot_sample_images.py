"""
=============
Sample Images
=============

This example shows off all of the sample images provided by SunPy.
"""
import matplotlib.pyplot as plt

from sunpy.data.sample import file_dict
import sunpy.map

###############################################################################
# We create a loop through all of the sample shortcuts and pick out the images
# that can be opened by Map. We then plot each one.
for key in file_dict:
    if key.count('IMAGE'):
        try:
            print(key)
            try:
                m = sunpy.map.Map(file_dict[key])
                plt.figure()
                m.plot()
                plt.show()
            except:
                pass
        except:
            pass
