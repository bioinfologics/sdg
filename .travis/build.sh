#! /bin/bash
# Bash Script that builds project

export PATH="${HOME}"/swig/bin:"${HOME}"/bin:"${HOME}"/cmake/bin:"${PATH}":"${HOME}"/Python/3.7/bin

mkdir build
cd build
cmake .. ${CMAKE_OPTIONS}
make all -j8
echo "" > ./docs/html/.nojekyll
echo "" > ./doc/sphinx/.nojekyll
make test

tar cz sdg-* libsdg.* pysdg/_pysdg.so pysdg/__init__.py pysdg/pysdg.py > sdg-${TRAVIS_OS_NAME}.tar.gz