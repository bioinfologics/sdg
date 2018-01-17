#! /bin/bash
# Bash Script that builds project
mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8

make test
ctest

(touch ./docs/html/.nojekyll || true)
