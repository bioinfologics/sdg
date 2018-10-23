#! /bin/bash
# Bash Script that builds project

export PATH="${HOME}"/swig/bin:"${PATH}"

mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8
echo "" > ./docs/html/.nojekyll
make test
