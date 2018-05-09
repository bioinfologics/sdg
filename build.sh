#! /bin/bash
# Bash Script that builds project
mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8
echo "" > ./docs/html/.nojekyll
make test
