#!/bin/bash
# Skip RTD build if commit message contains [skip ci] or [ci skip]
if git log -1 --format=%B | grep -iqE '\[(ci skip|skip ci)\]'; then
  exit 183
fi