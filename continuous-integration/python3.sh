#!/bin/bash
E=0

# This script runs a passing subset of tests under Python 3.

python setup.py install

python -c "import sunpy.data"
python -c "import sunpy.data; sunpy.data.download_sample_data()"
python -c "import sunpy.data.sample"

python setup.py test -P util
python setup.py test -P time
python setup.py test -P map
python setup.py test -P io
python setup.py test -P image
python setup.py test -P sun
python setup.py install || E=$?||$E
python -c "import sunpy.data" || E=$?||$E
python -c "import sunpy.data; sunpy.data.download_sample_data()" || E=$?||$E
python -c "import sunpy.data.sample" || E=$?||$E

py.test sunpy/time || E=$?||$E
py.test sunpy/map || E=$?||$E
py.test sunpy/io || E=$?||$E
py.test sunpy/image || E=$?||$E
py.test sunpy/sun || E=$?||$E

exit $E

