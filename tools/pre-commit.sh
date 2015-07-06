#!/bin/bash

PY_TEST=py.test
SPHINX=sphinx-build
BASE_DIR="$(git rev-parse --show-toplevel)"

NUM_SUNPY=$(git diff --name-only --stat HEAD $BASE_DIR/sunpy | wc -l)
NUM_DOCS=$(git diff --name-only --stat HEAD $BASE_DIR/doc | wc -l)
#IF code has changed run the tests:
if [ $NUM_SUNPY -gt  0 ]; then
    $PY_TEST $BASE_DIR
fi

#If the code OR the docs have changed run the docs
if  [ $NUM_SUNPY -gt 0 -o  $NUM_DOCS -gt  0 ]; then

    CWD="$(pwd)"
    DOC_DIR=$BASE_DIR/doc/source
    rm -r $DOC_DIR/../build
    rm -r $DOC_DIR/api
    rm -r $DOC_DIR/_generated

    cd $DOC_DIR
    $SPHINX -W -b html -d ../build/doctrees . ../build/html
    cd $CWD
fi
