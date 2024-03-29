cmake_minimum_required (VERSION 3.13)

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
if (POLICY CMP0028)
    cmake_policy(SET CMP0028 NEW)
endif (POLICY CMP0028)

if (POLICY CMP0078)
    cmake_policy(SET CMP0078 NEW)
endif (POLICY CMP0078)

if (POLICY CMP0086)
    cmake_policy(SET CMP0086 OLD)
endif (POLICY CMP0086)

if (APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)


PROJECT (sdg VERSION 0.1)

if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
    message(FATAL_ERROR "Unsupported compiler version, please use GCC version 6 or newer")
endif()

SET(CMAKE_RULE_MESSAGES OFF)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
ENABLE_TESTING()
OPTION(BUILD_TESTING "Build tests." OFF)
OPTION(BUILD_DOC "Build documentation" OFF)
OPTION(BUILD_SDGPYTHON "Build pybind11 SDGpython" ON)
OPTION(BUILD_COV "Build coverage" OFF)
SET(MAX_WORKSPACE_VECTOR_SIZE 100 CACHE STRING "Maximum number of items in a workspace vector")

set(CMAKE_CXX_STANDARD 14)

# Configure CCache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Fill the version variables in case the compilation isn't done from a git repo
if("${GIT_ORIGIN_URL}" STREQUAL "")
    # Get the remote url
    execute_process(
            COMMAND git remote get-url origin
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_ORIGIN_URL
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if("${GIT_ORIGIN_URL}" STREQUAL "")
        add_definitions("-DGIT_ORIGIN_URL=\"not set\"")
    else()
        add_definitions("-DGIT_ORIGIN_URL=\"${GIT_ORIGIN_URL}\"")
    endif()
else()
    add_definitions("-DGIT_ORIGIN_URL=\"${GIT_ORIGIN_URL}\"")
endif()

if("${GIT_COMMIT_HASH}" STREQUAL "")
    # Get the latest abbreviated commit hash of the working branch
    execute_process(
            COMMAND git log -1 --format=%h
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if("${GIT_COMMIT_HASH}" STREQUAL "")
        add_definitions("-DGIT_COMMIT_HASH=\"experimental\"")
    else()
        execute_process(
                COMMAND git log -1 --format=%h
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if("${GIT_STATUS}" STREQUAL "")
            # Check if the git directory is clean
            execute_process(
                    COMMAND git status -s
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_STATUS
                    OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            if("${GIT_STATUS}" STREQUAL "")
                add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
            else()
                add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH} (dirty)\"")
            endif()
        endif()
    endif()
    MESSAGE("STATUS - ${GIT_COMMIT_HASH}")
endif()


if("${GIT_BRANCH}" STREQUAL "")
    # Get the current working branch
    execute_process(
            COMMAND git rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if("${GIT_BRANCH}" STREQUAL "")
        add_definitions("-DGIT_BRANCH=\"not_from_git\"")
    else()
        add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")
    endif()
else()
    add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")
endif()

add_definitions("-DMAX_WORKSPACE_VECTOR_SIZE=${MAX_WORKSPACE_VECTOR_SIZE}")



## Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -Wno-deprecated")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fno-omit-frame-pointer -ftemplate-backtrace-limit=0 -Wall")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO"${CMAKE_CXX_FLAGS_DEBUG} -O3 -g -fno-omit-frame-pointer")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package (Threads)

include(FindOpenMP)
if(OPENMP_FOUND)
    #message("OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    #message("OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #message("OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
    set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${ZLIB_LIBRARIES})
endif()


## Include directories
include_directories("src/")
include_directories("deps/")
## Set sources variables

SET(header_files
    src/sdglib/graph/SequenceDistanceGraph.hpp
    src/sdglib/graph/SequenceDistanceGraphPath.hpp
    src/sdglib/graph/DistanceGraph.hpp
    src/sdglib/utilities/OutputLog.hpp
    src/sdglib/datastores/PairedReadsDatastore.hpp
    src/sdglib/datastores/LinkedReadsDatastore.hpp
    src/sdglib/datastores/LongReadsDatastore.hpp
    src/sdglib/datastores/KmerCounter.hpp
    src/sdglib/datastores/ReadSequenceBuffer.hpp
    src/sdglib/mappers/LongReadsMapper.hpp
    src/sdglib/mappers/LinkedReadsMapper.hpp
    src/sdglib/mappers/PairedReadsMapper.hpp
    src/sdglib/processors/LinkageUntangler.hpp
    src/sdglib/workspace/WorkSpace.hpp
    src/sdglib/processors/GraphEditor.hpp
    src/sdglib/processors/GraphMaker.hpp
    src/sdglib/indexers/UniqueKmerIndex.hpp
    src/sdglib/indexers/SatKmerIndex.hpp
    src/sdglib/indexers/NKmerIndex.hpp
    src/sdglib/views/NodeView.hpp
    src/sdglib/workspace/Journal.hpp
    src/sdglib/batch_counter/BatchKmersCounter.hpp
    deps/xxhash/xxhash.c
    src/sdglib/mappers/SequenceMapper.hpp
    src/sdglib/processors/LinkageMaker.hpp)

SET(source_files ${source_files}
    src/sdglib/graph/SequenceDistanceGraph.cc
    src/sdglib/graph/SequenceDistanceGraphPath.cc
    src/sdglib/graph/DistanceGraph.cc
    src/sdglib/utilities/OutputLog.cc
    src/sdglib/datastores/PairedReadsDatastore.cc
    src/sdglib/datastores/LinkedReadsDatastore.cc
    src/sdglib/datastores/LongReadsDatastore.cc
    src/sdglib/datastores/KmerCounter.cc
    src/sdglib/datastores/ReadSequenceBuffer.cc
    src/sdglib/mappers/LongReadsMapper.cc
    src/sdglib/mappers/LinkedReadsMapper.cc
    src/sdglib/mappers/PairedReadsMapper.cc
    src/sdglib/processors/LinkageUntangler.cc
    src/sdglib/workspace/WorkSpace.cc
    src/sdglib/processors/GraphEditor.cc
    src/sdglib/processors/GraphMaker.cc
    src/sdglib/indexers/UniqueKmerIndex.cc
    src/sdglib/indexers/SatKmerIndex.cc
    src/sdglib/indexers/NKmerIndex.cc
    src/sdglib/views/NodeView.cc
    src/sdglib/workspace/Journal.cc
    src/sdglib/batch_counter/BatchKmersCounter.cc
    deps/xxhash/xxhash.c
    src/sdglib/mappers/SequenceMapper.cc
    src/sdglib/processors/LinkageMaker.cc
    src/sdglib/processors/GraphContigger.cc
    src/sdglib/processors/GraphContigger.hpp
    src/sdglib/mappers/LongReadsRecruiter.cc
    src/sdglib/mappers/LongReadsRecruiter.hpp
    src/sdglib/processors/PathFinder.cc
    src/sdglib/processors/PathFinder.hpp
    src/sdglib/mappers/PerfectMatcher.cc
    src/sdglib/mappers/PerfectMatcher.hpp
    src/sdglib/processors/Strider.cc
    src/sdglib/processors/Strider.hpp
    src/sdglib/processors/GraphPatcher.cc
    src/sdglib/processors/GraphPatcher.hpp
    src/sdglib/views/TangleView.cc
    src/sdglib/views/TangleView.hpp
    src/sdglib/processors/LineFiller.h
    src/sdglib/processors/LineFiller.cc
    src/sdglib/processors/CountFilter.cc
    src/sdglib/processors/CountFilter.hpp
    src/sdglib/graph/ReadThreadsGraph.cc
    src/sdglib/graph/ReadThreadsGraph.hpp
    src/sdglib/processors/HappySorter.cc
    src/sdglib/processors/HappySorter.hpp
    src/sdglib/processors/TotalSorter.cc
    src/sdglib/processors/TotalSorter.hpp
    src/sdglib/processors/RTGCluster.cc
    src/sdglib/processors/RTGCluster.hpp
        src/sdglib/processors/RTGPartition.cc
        src/sdglib/processors/RTGPartition.hpp)

## Libraries
add_library(sdg SHARED
            ${source_files}
            )

target_link_libraries(sdg ${DEPENDENCY_LIBS})

add_executable(sdg-datastore
        src/sdg-datastore.cc
        )

target_link_libraries(sdg-datastore sdg)

MESSAGE("Building simple SDGpython")
FIND_PACKAGE(PythonInterp 3 REQUIRED)
FIND_PACKAGE(PythonLibs 3 REQUIRED)
ADD_LIBRARY( SDGpython SHARED interfaces/SDGpython.cc)
TARGET_INCLUDE_DIRECTORIES( SDGpython PUBLIC
        ${PYTHON_INCLUDE_DIRS}
        "interfaces/simple_python"
)
TARGET_LINK_LIBRARIES( SDGpython sdg )
# The library must not have any prefix and should be located in
# a subfolder that includes the package name. The setup will be
# more complicated otherwise.
SET_TARGET_PROPERTIES( SDGpython
        PROPERTIES
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/SDGpython"
        )

# This is required for linking the library under Mac OS X. Moreover,
# the suffix ensures that the module can be found by the interpreter
# later on.
IF( APPLE )
    SET_TARGET_PROPERTIES( SDGpython
            PROPERTIES
            LINK_FLAGS "-undefined dynamic_lookup"
            SUFFIX     ".so"
            )
ENDIF()

install(TARGETS
        sdg-datastore
        sdg
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(DIRECTORY src/sdglib DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")


