cmake_minimum_required (VERSION 2.8)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)
if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif (POLICY CMP0074)


if (APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)


PROJECT (bsg VERSION 0.1)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
ENABLE_TESTING()
OPTION(BUILD_TESTING "Build tests." ON)
OPTION(BUILD_DOC "Build documentation" OFF)
OPTION(BUILD_PYTHON_INTERFACE "Build Python interface" OFF)
OPTION(BUILD_R_INTERFACE "Build R interface" OFF)
OPTION(BUILD_COV "Build coverage" OFF)

set(CMAKE_CXX_STANDARD 14)

# Configure CCache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if (BUILD_PYTHON_INTERFACE)
    MESSAGE("Building docs for python interface")
    SET(BUILD_DOC ON)
endif ()
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
    add_definitions("-DGIT_COMMIT_HASH=\"not_from_git\"")
else()
    add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")
endif()

## Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -w -Wno-deprecated -std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -ftemplate-backtrace-limit=0 -Wall")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO"${CMAKE_CXX_FLAGS_DEBUG} -O3 -g -fno-omit-frame-pointer")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package (Threads)

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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

## Include directories
include_directories("src/")
include_directories("deps/cxxopts/")
include_directories("deps/kseqcpp/")
## Set sources variables

## Libraries
add_library(bsg SHARED
            src/sglib/graph/SequenceGraph.cc
            src/sglib/graph/SequenceSubGraph.cpp
            src/sglib/graph/SequenceGraphPath.cpp
            src/sglib/graph/LinkageDiGraph.cc
            src/sglib/filesystem/helpers.cc
            src/sglib/logger/OutputLog.cc
            src/sglib/datastores/PairedReadsDatastore.cc
            src/sglib/datastores/LinkedReadsDatastore.cc
            src/sglib/datastores/PathsDatastore.cc
            src/sglib/datastores/LongReadsDatastore.cpp
            src/sglib/mappers/LongReadMapper.cpp
            src/sglib/mappers/threader/NodeMapper.cpp
            src/sglib/mappers/threader/MappingThreader.cpp
            src/sglib/mappers/LinkedReadMapper.cc
            src/sglib/mappers/PairedReadMapper.cc
            src/sglib/processors/Untangler.cc
            src/sglib/processors/PathExplorer.cpp
            src/sglib/processors/LinkageUntangler.cc
            src/sglib/workspace/WorkSpace.cc
            src/sglib/processors/KmerCompressionIndex.cc
            src/sglib/processors/GraphEditor.cc
            src/sglib/processors/GraphMaker.cc
            src/sglib/processors/LocalHaplotypeAssembler.cc
            src/sglib/hash/xxhash.c
        )

add_library(bsg_static STATIC
        src/sglib/graph/SequenceGraph.cc
        src/sglib/graph/SequenceSubGraph.cpp
        src/sglib/graph/SequenceGraphPath.cpp
        src/sglib/graph/LinkageDiGraph.cc
        src/sglib/filesystem/helpers.cc
        src/sglib/logger/OutputLog.cc
        src/sglib/datastores/PairedReadsDatastore.cc
        src/sglib/datastores/LinkedReadsDatastore.cc
        src/sglib/datastores/PathsDatastore.cc
        src/sglib/datastores/LongReadsDatastore.cpp
        src/sglib/mappers/LongReadMapper.cpp
        src/sglib/mappers/threader/NodeMapper.cpp
        src/sglib/mappers/threader/MappingThreader.cpp
        src/sglib/mappers/LinkedReadMapper.cc
        src/sglib/mappers/PairedReadMapper.cc
        src/sglib/processors/Untangler.cc
        src/sglib/processors/PathExplorer.cpp
        src/sglib/processors/LinkageUntangler.cc
        src/sglib/workspace/WorkSpace.cc
        src/sglib/processors/KmerCompressionIndex.cc
        src/sglib/processors/GraphEditor.cc
        src/sglib/processors/GraphMaker.cc
        src/sglib/processors/LocalHaplotypeAssembler.cc
        src/sglib/hash/xxhash.c
        )

target_link_libraries(bsg ${DEPENDENCY_LIBS})
target_link_libraries(bsg_static ${DEPENDENCY_LIBS})

if (BUILD_COV)
    set_target_properties(bsg PROPERTIES
            COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
            LINK_FLAGS    "--coverage")
endif()

add_executable(bsg-datastore
        src/bsg-datastore.cc
        )

add_executable(bsg-kmerspectra
        src/bsg-kmerspectra.cc
        )

add_executable(bsg-workspace
        src/bsg-workspace.cc
        )

add_executable(bsg-mapper
        src/bsg-mapper.cc
        )

add_executable(bsg-untangler
        src/bsg-untangler.cc
        )

add_executable(bsg-lhapassembler
        src/bsg-lhapassembler.cc
        )

add_executable(seq-thread
        src/seq-thread.cpp
        )

add_executable(gfa-extract
        src/gfa-extract.cc)

add_executable(bsg-dbg
        src/bsg-dbg.cc)

target_link_libraries(bsg-workspace bsg_static)
target_link_libraries(bsg-datastore bsg_static)
target_link_libraries(bsg-kmerspectra bsg_static)
target_link_libraries(bsg-mapper bsg_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-untangler bsg_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg-lhapassembler bsg_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(gfa-extract bsg_static)
target_link_libraries(seq-thread bsg_static)
target_link_libraries(bsg-dbg bsg_static)


IF(BUILD_TESTING)
    include_directories("tests/include")
    add_library( CatchMain OBJECT tests/sglib/catch-main.cpp tests/include/catch.hpp )
    add_executable(runBasicTests
            $<TARGET_OBJECTS:CatchMain>
            tests/sglib/files_tests.cpp
            tests/sglib/smr_tests.cpp
            tests/sglib/kmerfactory_tests.cpp
            tests/sglib/index_tests.cpp
            )
    IF(BUILD_COV)
        set_target_properties(runBasicTests PROPERTIES
                COMPILE_FLAGS "-g -O0 -ftest-coverage -fprofile-arcs"
                LINK_FLAGS    "--coverage -lgcov")
    ENDIF()
    target_link_libraries(runBasicTests bsg)
    add_test(basicTests runBasicTests)
ENDIF()

IF(BUILD_PYTHON_INTERFACE)
    MESSAGE("Building python interface")
    FIND_PACKAGE(SWIG REQUIRED)
    INCLUDE(${SWIG_USE_FILE})
    set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/pybsg)
    configure_file(interfaces/python/__init__.py.in pybsg/__init__.py)
    configure_file(interfaces/SequenceGraph.i.in SequenceGraph.i)

    FIND_PACKAGE(PythonInterp 3 REQUIRED)
    FIND_PACKAGE(PythonLibs 3 REQUIRED)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY SWIG_MODULE_NAME pybsg)

    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTIES CPLUSPLUS ON)
    if (NOT APPLE)
        set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY COMPILE_DEFINITIONS "SWIGWORDSIZE64")
        set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY GENERATED_COMPILE_DEFINITIONS "SWIGWORDSIZE64")
    endif()

    if(DOXYGEN_FOUND)
        add_custom_command (
                OUTPUT ${CMAKE_BINARY_DIR}/python_docs.i
                COMMAND ${CMAKE_SOURCE_DIR}/doxy2swig.py
                --quiet ${CMAKE_BINARY_DIR}/docs/xml/index.xml
                ${CMAKE_BINARY_DIR}/python_docs.i
                COMMENT "Generating python_docs.i from Doxygen output for use by SWIG"
                DEPENDS doc_doxygen
                VERBATIM
        )
        add_custom_target(doxy2swig)
        SET_PROPERTY(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY SWIG_COMPILE_OPTIONS -O)
        SET_PROPERTY(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY DEPENDS ${CMAKE_BINARY_DIR}/python_docs.i)
    endif(DOXYGEN_FOUND)
    SWIG_ADD_LIBRARY(pybsg LANGUAGE python SOURCES ${CMAKE_BINARY_DIR}/SequenceGraph.i)
    SWIG_LINK_LIBRARIES(pybsg ${PYTHON_LIBRARIES} bsg_static)

    install(TARGETS
            _pybsg
            RUNTIME DESTINATION pybsg
            ARCHIVE DESTINATION pybsg
            LIBRARY DESTINATION pybsg)

ENDIF()

IF(BUILD_R_INTERFACE)
    MESSAGE("Building R interface")
    FIND_PACKAGE(SWIG REQUIRED)
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
    FIND_PACKAGE(R REQUIRED)
    INCLUDE_DIRECTORIES(${R_INCLUDE_DIR})
    INCLUDE(${SWIG_USE_FILE})
    set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/rbsg)
    configure_file(interfaces/SequenceGraph.i.in SequenceGraph.i)

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTIES CPLUSPLUS ON)
    set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceGraph.i PROPERTY SWIG_MODULE_NAME rbsg)
    SWIG_ADD_LIBRARY(rbsg LANGUAGE r SOURCES ${CMAKE_BINARY_DIR}/SequenceGraph.i)
    SWIG_LINK_LIBRARIES(rbsg ${R_LIBRARIES} bsg_static)

    install(TARGETS
            rbsg
            RUNTIME DESTINATION rbsg
            ARCHIVE DESTINATION rbsg
            LIBRARY DESTINATION rbsg)
ENDIF()

install(TARGETS
        bsg-workspace
        bsg-datastore
        bsg-kmerspectra
        bsg-mapper
        bsg-untangler
        bsg-lhapassembler
        gfa-extract
        seq-thread
        bsg
        bsg_static
        runBasicTests
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(DIRECTORY src/sglib
        DESTINATION include)
