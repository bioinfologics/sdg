cmake_minimum_required (VERSION 3.0)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

if (APPLE)
       set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)


PROJECT (bsg)
ENABLE_TESTING()
OPTION(BUILD_TESTING "Build tests." ON)
OPTION(BUILD_DOC "Build documentation" OFF)
OPTION(BUILD_COV "Build coverage" OFF)

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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")


find_package (Threads)

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#add_definitions(${CMAKE_CXX_FLAGS})

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
    set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${ZLIB_LIBRARIES})
endif()

find_package(BZip2)
if (BZIP2_FOUND)
    add_definitions(-DHAVE_BZIP2)
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${BZIP2_LIBRARIES})
endif()

# check if Doxygen is installed
IF(BUILD_DOC)
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
ENDIF()

if(BUILD_TESTING)
    add_subdirectory("tests/googletest")
    include_directories("tests/googletest/googletest/include/")
    include_directories("tests/googletest/googletest/src/")
ENDIF()
## Include directories
include_directories("src/")
include_directories("deps/cxxopts/include")
include_directories("deps/kseqcpp/")
## Set sources variables

## Libraries
add_library(sequence_graph SHARED
            src/sglib/datastores/LinkedReadsDatastore.cc
            src/sglib/factories/StrandedMinimiserSketchFactory.cc
            src/sglib/filesystem/check_or_create_directory.cc
            src/sglib/logger/OutputLog.cc
            src/sglib/mappers/LinkedReadMapper.cc
            src/sglib/processors/TagWalker.cc
            src/sglib/processors/Untangler.cc
            src/sglib/GraphPartitioner.cc
            src/sglib/HaplotypeScorer.cpp
            src/sglib/KmerCompressionIndex.cc
            src/sglib/PairedReadMapper.cc
            src/sglib/PhaseScaffolder.cpp
            src/sglib/Scaffolder.cc
            src/sglib/SequenceGraph.cc
            src/sglib/SequenceGraphJournal.cc
            src/sglib/SequenceThreader.cpp
            src/sglib/UniqueKmerIndex.cpp
            src/sglib/WorkSpace.cc
            src/sglib/processors/PathExplorer.cpp
            src/sglib/aligners/algorithms/ksw2_extz2_sse.c
            )

target_link_libraries(sequence_graph ${DEPENDENCY_LIBS})

if (BUILD_COV)
    set_target_properties(sequence_graph PROPERTIES
                      COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
                      LINK_FLAGS    "--coverage")
endif()

target_link_libraries(sequence_graph ${DEPENDENCY_LIBS})

install(TARGETS sequence_graph DESTINATION lib)
install(FILES
            src/utilities/matrix.h
            src/sglib/datastores/LinkedReadsDatastore.hpp
            src/sglib/factories/ContigBlockFactory.h
            src/sglib/factories/ContigLink.h
            src/sglib/factories/KMerFactory.h
            src/sglib/factories/KMerIDXFactory.h
            src/sglib/factories/StrandedMinimiserSketchFactory.h
            src/sglib/filesystem/check_or_create_directory.h
            src/sglib/logger/OutputLog.h
            src/sglib/mappers/LinkedReadMapper.hpp
            src/sglib/processors/TagWalker.hpp
            src/sglib/processors/Untangler.hpp
            src/sglib/readers/Common.h
            src/sglib/readers/FileReader.h
            src/sglib/readers/SequenceGraphReader.h
            src/sglib/GraphPartitioner.hpp
            src/sglib/HaplotypeScorer.hpp
            src/sglib/KmerCompressionIndex.hpp
            src/sglib/PairedReadMapper.hpp
            src/sglib/PhaseScaffolder.h
            src/sglib/Scaffolder.hpp
            src/sglib/SequenceGraph.hpp
            src/sglib/SequenceGraphJournal.hpp
            src/sglib/SequenceThreader.h
            src/sglib/SMR.h
            src/sglib/UniqueKmerIndex.h
            src/sglib/WorkSpace.hpp
            src/sglib/aligners/algorithms/SmithWaterman.h
            src/sglib/aligners/submat/SubstitutionMatrix.h
            src/sglib/processors/PathExplorer.h
            src/sglib/aligners/algorithms/ksw2.h
        DESTINATION include
        )



add_executable(bsg-datastore
        src/bsg-datastore.cc
        )

add_executable(bsg-kmerspectra
        src/bsg-kmerspectra.cc
        )

add_executable(bsg-workspace
        src/bsg-workspace.cc
        )

add_executable(bsg-untangler
        src/bsg-untangler.cc
        )

## add executables and respective sources
add_executable(gfa-qc
        src/gfa-qc.cc)

## add executables and respective sources
add_executable(gfa-dumppath
        src/gfa-dumppath.cc)

## add executables and respective sources
add_executable(gfa-pairscaff
        src/gfa-pairscaff.cc)

add_executable(gfa-bfphaser
        src/gfa-bfphaser.cc)

add_executable(phaser
               src/phaser.cpp)

add_executable(gfa-align
               src/gfa-align.cc)

add_executable(kc_comp
               src/comp_kc.cc)

add_executable(seq-sorter
               src/seq-sorter.cpp)

add_executable(seq-thread
        src/seq-thread.cpp)

add_executable(map-lr
               src/map-lr.cc)

add_executable(gfa-extract
        src/gfa-extract.cc)

target_link_libraries(bsg-workspace sequence_graph)
target_link_libraries(bsg-datastore sequence_graph)
target_link_libraries(gfa-dumppath sequence_graph)
target_link_libraries(gfa-qc sequence_graph)
target_link_libraries(seq-sorter sequence_graph)
target_link_libraries(gfa-pairscaff sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-bfphaser sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-untangler sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-align sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(phaser sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(seq-thread sequence_graph)
target_link_libraries(phaser sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(map-lr sequence_graph ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-extract sequence_graph)

IF(BUILD_TESTING)
    add_executable(runBasicTests tests/sglib/basic_check.cc)
    IF(BUILD_COV)
        set_target_properties(runBasicTests PROPERTIES
                          COMPILE_FLAGS "-g -O0 -ftest-coverage -fprofile-arcs"
                          LINK_FLAGS    "--coverage")
    ENDIF()
    target_link_libraries(runBasicTests gtest gtest_main sequence_graph)
    add_test(basicTests runBasicTests)
ENDIF()

install(TARGETS
        gfa-dumppath
        gfa-qc
        seq-sorter
        gfa-pairscaff
        gfa-align
        phaser
        seq-thread
        runBasicTests
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib/static)
install(DIRECTORY src/sglib
        DESTINATION include)
