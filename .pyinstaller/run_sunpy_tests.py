import os
import sys

sys.exit(os.system('pytest -p no:warnings --doctest-rst -m "not mpl_image_compare" --pyargs sunpy'))
