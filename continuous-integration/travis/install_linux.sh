###############################################################################
# Install Linux Packages
###############################################################################

# Installation of non-Python dependencies that are in the travis approved
# packages list is now in .travis.yml

# Install more upto date openjpeg library.
wget http://openjpeg.googlecode.com/files/openjpeg-1.5.0-Linux-x86_64.tar.gz
tar xvzf openjpeg-1.5.0-Linux-x86_64.tar.gz --strip-components=1
export PATH=$HOME:$PATH

export LIBRARY_PATH=$(pwd)/lib
export LD_LIBRARY_PATH=$(pwd)/lib

###############################################################################
# Install miniconda
###############################################################################
MINICONDA_URL="http://repo.continuum.io/miniconda"
MINICONDA_FILE="Miniconda-3.5.5-Linux-x86_64.sh"
wget "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

export PATH=$HOME/miniconda/bin:$PATH

conda update --yes conda
conda install -q --yes binstar conda-build
