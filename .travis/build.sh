#! /bin/bash
# Bash Script that builds project

export PATH="${HOME}"/swig/bin:"${HOME}"/bin:"${HOME}"/cmake/bin:"${PATH}"

mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8
echo "" > ./docs/html/.nojekyll
make test
