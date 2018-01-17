#! /bin/bash
# Bash Script that builds project
mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make -j8

find . -iname "*.gcda"

PATH=$PATH:./bin make test
PATH=$PATH:./bin ctest

find . -iname "*.gcda"

touch ./docs/html/.nojekyll || true