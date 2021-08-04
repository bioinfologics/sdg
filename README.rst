|license|

.. |license| image:: https://img.shields.io/badge/license-MIT-green.svg
    :target: https://github.com/bioinfologics/bsg/blob/master/LICENSE

Sequence Distance Graph
========================

The Sequence Distance Graph (**SDG**) is a framework to work with genome graphs and sequencing data. It provides a workspace built around a Sequence Distance Graph, datastores for paired, linked and long reads, read mappers, and k-mer counters. It can be used to perform different types of sequence analyses.

SDG can be run on Linux and MacOS, and requires enough RAM to hold the WorkSpace completely in memory, which will depend on the dataset. Space to hold the uncompressed sequences on the datastores on disk will also be required.

SDG is used as a Python module through its pybind interface.

SDG has been published on an F1000 software article: https://f1000research.com/articles/8-1490/v1.

Installing SDG releases
#######################

**Notice: this refers to older release version of SDG, we are preparing the new version for release and will update this page soon. In the meantime, please use these instructions with the release code.**

SDG can be installed via pre-compiled binaries from https://github.com/bioinfologics/sdg/releases. The binaries have been built using Python3 and GCC version 6 from the Ubuntu package manager for the Linux version. The MacOS version dependencies were obtained using Homebrew (Python3, GCC-6 and SWIG).

Compiling the release versions requires the following dependencies:

- GCC 6+
- CMake 3.14.3+
- Python 3.7+ (For Python interface)
- Git LFS (For tests)

The installation process consists of generating the configuration files using CMake and using make. By default make install will install to system wide directories, to specify a different directory use `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/v3.13/variable/CMAKE_INSTALL_PREFIX.html#cmake-install-prefix>`_.

.. code-block:: bash

    wget https://github.com/bioinfologics/sdg/archive/v1.0_rc8.tar.gz
    cd sdg-1.0_rc8
    mkdir build
    cd build
    cmake -DBUILD_SIMPLE_PYTHON_INTERFACE=on ../
    make
    make test
    make install


Usage
#####

Please refer to the publication: https://f1000research.com/articles/8-1490/v1.

Citing SDG
#######
Yanes L, Garcia Accinelli G, Wright J et al. A Sequence Distance Graph framework for genome assembly and analysis [version 1; peer review: 2 approved]. F1000Research 2019, 8:1490
(https://doi.org/10.12688/f1000research.20233.1)