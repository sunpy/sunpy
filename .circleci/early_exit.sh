#!/bin/bash

# Get the most recent commit message
commit_message=$(git log --pretty=%B -n 1)

# Check for skip directives in the commit message
if [[ "$commit_message" == *"[ci skip]"* || "$commit_message" == *"[skip ci]"* ]]; then
    echo "Skipping build because '[ci skip]' or '[skip ci]' was found in the commit message."
    circleci step halt
fi
