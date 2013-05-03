#!/bin/bash
#
# This is a simple script that I built to add, commit and re-install the sunpy library
# Author: Michael Malocha
# Date: 5-2-2013
#

NOW=$(date +"%F")
git add .
git commit -m "Test $NOW"
sudo python setup.py install