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

project (sg)

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

## Include directories
include_directories("src/")

## Set sources variables

## add executables and respective sources
add_executable(gfa-qc
        src/gfa-qc.cc
        src/sglib/SequenceGraph.cc src/sglib/SequenceGraph.hpp)

## add executables and respective sources
add_executable(gfa-dumppath
        src/gfa-dumppath.cc
        src/sglib/SequenceGraph.cc src/sglib/SequenceGraph.hpp)

## add executables and respective sources
add_executable(gfa-pairscaff
        src/gfa-pairscaff.cc
        src/sglib/FileReader.h
        src/sglib/KMerFactory.h
        src/sglib/KMerIDX.h
        src/sglib/SMR.h
        src/sglib/SequenceGraph.hpp
        src/sglib/SequenceGraph.cc
        src/sglib/PairedReadMapper.hpp
        src/sglib/PairedReadMapper.cc
        src/sglib/Scaffolder.cc src/sglib/Scaffolder.hpp)


target_link_libraries (gfa-pairscaff ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})
