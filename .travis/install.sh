#!/usr/bin/env bash
set -x
set -e
if [ "${TRAVIS_OS_NAME}" == linux ]; then
    # Need SWIG >= 3.0.8
    cd /tmp/ &&
    wget https://github.com/swig/swig/archive/rel-3.0.12.tar.gz &&
    tar zxf rel-3.0.12.tar.gz && cd swig-rel-3.0.12 &&
    ./autogen.sh && ./configure --prefix "${HOME}"/swig/ 1>/dev/null &&
    make >/dev/null &&
    make install >/dev/null;


    export DOXYGEN_VER=doxygen-1.8.14
    export DOXYGEN_TAR=${DOXYGEN_VER}.linux.bin.tar.gz
    export DOXYGEN_URL="http://ftp.stack.nl/pub/users/dimitri/${DOXYGEN_TAR}"
    wget -O - "${DOXYGEN_URL}" | tar xz -C ${HOME} bin/doxygen
fi
