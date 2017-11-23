cmake_minimum_required (VERSION 2.8)

IF (MY_COMPILER_OPTION STREQUAL "Intel")
    FIND_PROGRAM(ICPC_PROGRAM icpc ...)
    IF(ICPC_PROGRAM)
        MESSAGE(STATUS "** Compiling with Intel settings **")
        IF(ENV{CXX})
            MESSAGE(WARNING "Overwriting CXX envvar")
        ENDIF()
        SET(ENV{CXX} "${ICPC_PROGRAM}")
        SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
        SET(CMAKE_CXX_FLAGS_DEBUG "-g")
    ELSE()
        MESSAGE(FATAL_ERROR "Intel compiler not found")
    ENDIF()
ENDIF ()

if (APPLE)
       set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)


project (bsg)

## Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -w -Wno-deprecated -std=c++14")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


find_package (Threads)

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_definitions(${CMAKE_CXX_FLAGS})

## To move the bins to the final location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

find_package(BZip2 REQUIRED)
if (BZIP2_FOUND)
    add_definitions(-DHAVE_BZIP2)
    include_directories(${BZIP2_INCLUDE_DIRS})
endif()

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target(doc_doxygen ALL
                      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)


## Include directories
include_directories("src/")
include_directories("deps/cxxopts/include")
include_directories("deps/kseqcpp/")

## Set sources variables

## Libraries
add_library(sequence_graph SHARED
            src/sglib/SMR.h
            src/sglib/readers/Common.h
            src/sglib/readers/FileReader.h
            src/sglib/readers/SequenceGraphReader.h
            src/sglib/factories/ContigBlockFactory.h
            src/sglib/factories/ContigLink.h
            src/sglib/factories/KMerFactory.h
            src/sglib/factories/KMerIDXFactory.h
            src/sglib/SequenceGraph.hpp
            src/sglib/SequenceGraph.cc
            src/sglib/PairedReadMapper.hpp
            src/sglib/PairedReadMapper.cc
            src/sglib/HaplotypeScorer.hpp
            src/sglib/HaplotypeScorer.cpp
            src/sglib/SequenceGraphJournal.hpp
            src/sglib/SequenceGraphJournal.cc
            src/sglib/SequenceGraph.cc
            src/sglib/PairedReadMapper.cc
            src/sglib/Scaffolder.cc
            src/sglib/filesystem/check_or_create_directory.cc
            src/sglib/filesystem/check_or_create_directory.h
            src/sglib/PhaseScaffolder.cpp
            src/sglib/PhaseScaffolder.h)

target_link_libraries(sequence_graph ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})

## add executables and respective sources
add_executable(gfa-qc
        src/gfa-qc.cc)

## add executables and respective sources
add_executable(gfa-dumppath
        src/gfa-dumppath.cc)

## add executables and respective sources
add_executable(gfa-pairscaff
        src/gfa-pairscaff.cc)

add_executable(phaser
               src/phaser.cpp)

add_executable(gfa-align
               src/gfa-align.cc)

add_executable(kc_comp
               src/comp_kc.cc)

add_executable(seq-sorter
               src/seq-sorter.cpp)

target_link_libraries(gfa-dumppath sequence_graph)
target_link_libraries(gfa-qc sequence_graph)
target_link_libraries(seq-sorter sequence_graph)
target_link_libraries(gfa-pairscaff sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-align sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (phaser sequence_graph ${CMAKE_THREAD_LIBS_INIT})
