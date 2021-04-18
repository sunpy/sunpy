import os
import sys
import pytest
import shutil

sys.exit(os.system('pytest -p no:warnings --doctest-rst -m "not mpl_image_compare" --pyargs sunpy'))                              