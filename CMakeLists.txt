cmake_minimum_required (VERSION 2.8)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()


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

# Get the remote url
execute_process(
        COMMAND git remote get-url origin
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_ORIGIN_URL
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check if the git directory is clean
execute_process(
        COMMAND git status -s
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Fill the version variables in case the compilation isn't done from a git repo
if("${GIT_ORIGIN_URL}" STREQUAL "")
    add_definitions("-DGIT_ORIGIN_URL=\"not set\"")
else()
    add_definitions("-DGIT_ORIGIN_URL=\"${GIT_ORIGIN_URL}\"")
endif()
if("${GIT_COMMIT_HASH}" STREQUAL "")
    add_definitions("-DGIT_COMMIT_HASH=\"experimental\"")
else()
    if("${GIT_STATUS}" STREQUAL "")
        add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
    else()
        add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH} (dirty)\"")
    endif()
endif()
if("${GIT_BRANCH}" STREQUAL "")
    add_definitions("-DGIT_COMMIT_HASH=\"not from git\"")
else()
    add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")
endif()

## Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -w -Wno-deprecated -std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fno-omit-frame-pointer")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO"${CMAKE_CXX_FLAGS_DEBUG} -O3 -g -fno-omit-frame-pointer")


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

    add_custom_target(doc_doxygen
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
        src/sglib/filesystem/check_or_create_directory.cc
        src/sglib/filesystem/check_or_create_directory.h
        src/sglib/KmerCompressionIndex.cc
        src/sglib/KmerCompressionIndex.hpp
        src/sglib/logger/OutputLog.h
        src/sglib/logger/OutputLog.cc
        src/sglib/datastores/LinkedReadsDatastore.cc
        src/sglib/datastores/LinkedReadsDatastore.hpp
        src/sglib/datastores/PairedReadsDatastore.cc
        src/sglib/datastores/PairedReadsDatastore.hpp
        src/sglib/mappers/ReadMapping.hpp
        src/sglib/mappers/LinkedReadMapper.cc
        src/sglib/mappers/LinkedReadMapper.hpp
        src/sglib/mappers/PairedReadMapper.cc
        src/sglib/mappers/PairedReadMapper.hpp
        src/sglib/datastores/PathsDatastore.cc
        src/sglib/datastores/PathsDatastore.hpp
        src/sglib/WorkSpace.cc
        src/sglib/WorkSpace.hpp
        src/sglib/LinkageDiGraph.cc
        src/sglib/LinkageDiGraph.hpp)

target_link_libraries(sequence_graph ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})

add_executable(bsg-datastore
        src/bsg-datastore.cc
        src/sglib/datastores/LinkedReadsDatastore.hpp
        src/sglib/datastores/LinkedReadsDatastore.cc
        src/sglib/logger/OutputLog.h
        src/sglib/logger/OutputLog.cc)

add_executable(bsg-kmerspectra
        src/bsg-kmerspectra.cc
        src/sglib/SequenceGraph.hpp
        src/sglib/SequenceGraph.cc
        src/sglib/filesystem/check_or_create_directory.h
        src/sglib/filesystem/check_or_create_directory.cc
        src/sglib/KmerCompressionIndex.hpp
        src/sglib/KmerCompressionIndex.cc
        src/sglib/logger/OutputLog.h
        src/sglib/logger/OutputLog.cc
        )
target_link_libraries(bsg-kmerspectra ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})

add_executable(bsg-workspace
        src/bsg-workspace.cc
        )
target_link_libraries(bsg-kmerspectra ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})

add_executable(bsg-mapper
        src/bsg-mapper.cc
        )

add_executable(bsg-untangler
        src/bsg-untangler.cc
        src/sglib/processors/Untangler.hpp
        src/sglib/processors/Untangler.cc
        src/sglib/processors/TagWalker.cc src/sglib/processors/TagWalker.hpp
        src/sglib/processors/FlowFollower.cc src/sglib/processors/FlowFollower.hpp
        src/sglib/processors/LinkageUntangler.cc src/sglib/processors/LinkageUntangler.hpp
        src/sglib/processors/GraphEditor.cc src/sglib/processors/GraphEditor.hpp
        src/sglib/processors/GraphMaker.cc src/sglib/processors/GraphMaker.hpp
        src/sglib/processors/LocalHaplotypeAssembler.cc src/sglib/processors/LocalHaplotypeAssembler.hpp)

add_executable(bsg-dbg
        src/bsg-dbg.cc
        src/sglib/processors/GraphMaker.cc src/sglib/processors/GraphMaker.hpp
        src/sglib/processors/GraphEditor.cc src/sglib/processors/GraphEditor.hpp)


add_executable(bsg-flowmaker
        src/bsg-flowmaker.cc
        src/sglib/processors/Untangler.hpp
        src/sglib/processors/Untangler.cc
        src/sglib/processors/TagWalker.cc
        src/sglib/processors/TagWalker.hpp
        src/sglib/processors/FlowFollower.cc
        src/sglib/processors/FlowFollower.hpp)

add_executable(test_graphmaker EXCLUDE_FROM_ALL
        src/test_graphmaker.cc
        src/sglib/processors/GraphMaker.cc
        src/sglib/processors/GraphMaker.hpp)



## add executables and respective sources
add_executable(gfa-qc EXCLUDE_FROM_ALL
        src/gfa-qc.cc)

## add executables and respective sources
add_executable(gfa-dumppath EXCLUDE_FROM_ALL
        src/gfa-dumppath.cc)

## add executables and respective sources
add_executable(gfa-pairscaff EXCLUDE_FROM_ALL
        src/gfa-pairscaff.cc)

add_executable(gfa-align EXCLUDE_FROM_ALL
               src/gfa-align.cc)

add_executable(kc_comp EXCLUDE_FROM_ALL
               src/comp_kc.cc)

add_executable(seq-sorter EXCLUDE_FROM_ALL
               src/seq-sorter.cpp)

target_link_libraries(bsg-workspace sequence_graph)
target_link_libraries(bsg-datastore sequence_graph)
target_link_libraries(bsg-kmerspectra sequence_graph)
target_link_libraries(bsg-mapper sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-flowmaker sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-untangler sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-dbg sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(test_graphmaker sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-dumppath sequence_graph)
target_link_libraries(gfa-qc sequence_graph)
target_link_libraries(seq-sorter sequence_graph)
target_link_libraries(gfa-pairscaff sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-align sequence_graph ${CMAKE_THREAD_LIBS_INIT})
