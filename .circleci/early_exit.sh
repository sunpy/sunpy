#!/bin/bash

commitmessage=$(git log --pretty=%B -n 1)
if [[ $commitmessage = *"[ci skip]"* ]]; then
    echo "Skipping build because [ci skip] found in commit message"
    circleci step halt
fi
