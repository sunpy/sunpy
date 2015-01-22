MINICONDA_URL="http://repo.continuum.io/miniconda"
MINICONDA_FILE="Miniconda3-3.7.3-MacOSX-x86_64.sh"
wget "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

export PATH=/Users/travis/miniconda3/bin:$PATH

conda update --yes conda
conda install -q --yes binstar conda-build
