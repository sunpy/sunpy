#!/bin/bash

NOW=$(date +"%F")
git add .
git commit -m "Test $NOW"
sudo python setup.py install