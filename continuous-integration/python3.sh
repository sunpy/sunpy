#!/bin/bash
E=0

# This script runs a passing subset of tests under Python 3.

python setup.py install || E=$?||$E
python -c "import sunpy.data" || E=$?||$E
python -c "import sunpy.data; sunpy.data.download_sample_data()" || E=$?||$E
python -c "import sunpy.data.sample" || E=$?||$E

py.test sunpy/time || E=$?||$E
py.test sunpy/map -m "not figure"|| E=$?||$E
py.test sunpy/io || E=$?||$E
py.test sunpy/image || E=$?||$E
py.test sunpy/sun || E=$?||$E
py.test sunpy/util || E=$?||$E
py.test sunpy/spectra || E=$?||$E
py.test sunpy/net || E=$?||$E

exit $E

