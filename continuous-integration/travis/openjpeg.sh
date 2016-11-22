#!/bin/bash

# Create the glymur config file to set the path for the openjpeg library
mkdir -p $HOME/.config/glymur
cat << EOF > $HOME/.config/glymur/glymurrc
[library]
openjp2:  $HOME/miniconda/envs/test/lib/libopenjp2.so
EOF
