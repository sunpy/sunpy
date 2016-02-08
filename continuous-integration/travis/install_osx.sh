# Guide at http://conda.pydata.org/docs/travis.html#the-travis-yml-file

MINICONDA_URL="http://repo.continuum.io/miniconda"
MINICONDA_FILE="Miniconda3-latest-MacOSX-x86_64.sh"
wget "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b -p $HOME/miniconda

export PATH=$HOME/miniconda/bin:$PATH

hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda install -q --yes binstar conda-build
