###############################################################################
# Install Linux Packages
###############################################################################
sudo apt-get upgrade
sudo apt-get install -qq libatlas-dev liblapack-dev gfortran
if [[ $TEST_MODE == 'sphinx' ]]; then sudo apt-get install graphviz texlive-latex-extra dvipng; fi

# Install more upto date openjpeg library.
wget http://openjpeg.googlecode.com/files/openjpeg-1.5.0-Linux-x86_64.tar.gz
sudo tar -xvf openjpeg-1.5.0-Linux-x86_64.tar.gz --strip-components=1 -C /

###############################################################################
# Install miniconda
###############################################################################
MINICONDA_URL="http://repo.continuum.io/miniconda"
MINICONDA_FILE="Miniconda-3.5.5-Linux-x86_64.sh"
wget "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

export PATH=$HOME/miniconda/bin:$PATH

conda update --yes conda
