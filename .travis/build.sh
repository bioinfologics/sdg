#! /bin/bash
# Bash Script that builds project

export PATH="${HOME}"/swig/bin:"${HOME}"/bin:"${HOME}"/cmake/bin:"${PATH}":"${HOME}"/Python/3.7/bin

mkdir build
cd build
mkdir product
cmake -DCMAKE_INSTALL_PREFIX=product .. ${CMAKE_OPTIONS}
make all -j8
echo "" > ./docs/html/.nojekyll
echo "" > ./doc/sphinx/.nojekyll
make test
make install

tar cz product > sdg-${TRAVIS_OS_NAME}.tar.gz