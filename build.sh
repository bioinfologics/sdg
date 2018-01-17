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
lcov --directory . --capture --output-file coverage.info
lcov --remove coverage.info '/usr/*' --output-file coverage.info
lcov --remove coverage.info '*Applications*' --output-file coverage.info
lcov --remove coverage.info '*googletest*' --output-file coverage.info
lcov --remove coverage.info '*v1*' --output-file coverage.info
lcov --list coverage.info
bash <(curl -s https://codecov.io/bash)