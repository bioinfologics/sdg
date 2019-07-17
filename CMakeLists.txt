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

SET(CMAKE_RULE_MESSAGES OFF)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
INCLUDE(FindPythonModule)
ENABLE_TESTING()
OPTION(BUILD_TESTING "Build tests." ON)
OPTION(BUILD_DOC "Build documentation" OFF)
OPTION(BUILD_PYTHON_INTERFACE "Build Python interface" OFF)
OPTION(BUILD_R_INTERFACE "Build R interface" OFF)
OPTION(BUILD_JULIA_INTERFACE "Build Julia interface" OFF)
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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fno-omit-frame-pointer -ftemplate-backtrace-limit=0 -Wall")

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


## Include directories
include_directories("src/")
include_directories("deps/cxxopts/")
include_directories("deps/kseqcpp/")
include_directories("deps/xxhash/")
## Set sources variables
SET(source_files ${source_files}
    src/sdglib/graph/SequenceDistanceGraph.cc
    src/sdglib/graph/SequenceSubGraph.cc
    src/sdglib/graph/SequenceGraphPath.cc
    src/sdglib/graph/DistanceGraph.cc
        src/sdglib/utilities/OutputLog.cc
    src/sdglib/datastores/PairedReadsDatastore.cc
    src/sdglib/datastores/LinkedReadsDatastore.cc
    src/sdglib/datastores/LongReadsDatastore.cc
    src/sdglib/datastores/KmerCounts.cc
    src/sdglib/datastores/ReadPathsDatastore.cc
    src/sdglib/datastores/ReadSequenceBuffer.cc
    src/sdglib/mappers/LongReadsMapper.cc
    src/sdglib/mappers/LinkedReadsMapper.cc
    src/sdglib/mappers/PairedReadsMapper.cc
        src/sdglib/processors/LinkageUntangler.cc
    src/sdglib/workspace/WorkSpace.cc
        src/sdglib/processors/GraphEditor.cc
    src/sdglib/processors/GraphMaker.cc
    src/sdglib/processors/LocalHaplotypeAssembler.cc
    src/sdglib/processors/HaplotypeConsensus.cc
    src/sdglib/indexers/UniqueKmerIndex.cc
    src/sdglib/indexers/SatKmerIndex.cc
    src/sdglib/indexers/NKmerIndex.cc
    src/sdglib/views/NodeView.cc
    src/sdglib/journal/OperationJournal.cc
    src/sdglib/batch_counter/BatchKmersCounter.cc
    deps/xxhash/xxhash.c
        src/sdglib/mappers/SequenceMapper.cc src/sdglib/mappers/SequenceMapper.hpp)
## Libraries
add_library(sdg SHARED
            ${source_files}
            )

add_library(sdg_static STATIC
            ${source_files}
            )

target_link_libraries(sdg ${DEPENDENCY_LIBS})
target_link_libraries(sdg_static ${DEPENDENCY_LIBS})

if (BUILD_COV)
    set_target_properties(sdg PROPERTIES
            COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
            LINK_FLAGS    "--coverage")
endif()

add_executable(sdg-datastore
        src/sdg-datastore.cc
        )

add_executable(sdg-kmercounts
               src/sdg-kmercounts.cc
        )

add_executable(sdg-workspace
        src/sdg-workspace.cc
        )

add_executable(sdg-mapper
        src/sdg-mapper.cc
        )

add_executable(sdg-dbg
        src/sdg-dbg.cc)

target_link_libraries(sdg-workspace sdg_static)
target_link_libraries(sdg-datastore sdg_static)
target_link_libraries(sdg-kmercounts sdg_static)
target_link_libraries(sdg-mapper sdg_static)
target_link_libraries(sdg-dbg sdg_static)


IF(BUILD_TESTING)
    enable_testing()
    include_directories("tests/include")
    add_library(CatchMain OBJECT tests/sdglib/catch-main.cc tests/include/catch.hpp)
    add_executable(runBasicTests
                   $<TARGET_OBJECTS:CatchMain>
                   tests/sdglib/files_tests.cc
                   tests/sdglib/kmerfactory_tests.cc
                   tests/sdglib/index_tests.cc
                   tests/sdglib/workspace_tests.cc tests/sdglib/batchCounter_tests.cc)
    IF(BUILD_COV)
        set_target_properties(runBasicTests PROPERTIES
                COMPILE_FLAGS "-g -O0 -ftest-coverage -fprofile-arcs"
                LINK_FLAGS    "--coverage -lgcov")
    ENDIF()
    target_link_libraries(runBasicTests sdg)
    add_test(NAME basicTests COMMAND $<TARGET_FILE:runBasicTests>)
ENDIF()

# check if Doxygen is installed
IF(BUILD_DOC)
    find_package(Doxygen REQUIRED)
    if (DOXYGEN_FOUND)
        get_target_property(SDG_PUBLIC_HEADER_DIR sdg INTERFACE_INCLUDE_DIRECTORIES)
        file(GLOB_RECURSE SDG_PUBLIC_HEADERS ${SDG_PUBLIC_HEADER_DIR}/*.h)
        # set input and output files
        set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})
        set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)
        set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")
        # note the option ALL which allows to build the docs together with the application
        add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                           DEPENDS ${SDG_PUBLIC_HEADERS} README.md src/sdglib/README.md src/sdglib/workspace/README.md
                           COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                           MAIN_DEPENDENCY Doxyfile
                           COMMENT "Generating docs")
        add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_OUT})

        FIND_PACKAGE(PythonInterp 3)
        FIND_PACKAGE(PythonLibs 3)

        find_package(Sphinx)
        if(SPHINX_FOUND AND NOT PYTHON_EXECUTABLE EQUAL "")
            find_python_module(breathe)
            if(PY_BREATHE)
                set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx)
                set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/doc/sphinx)
                set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

                add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                                   COMMAND
                                   ${SPHINX_EXECUTABLE} -b html
                                   # Tell Breathe where to find the Doxygen output
                                   -Dbreathe_projects.SDG=${DOXYGEN_OUTPUT_DIR}/xml
                                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                                   DEPENDS
                                   # Other docs files you want to track should go here (or in some variable)
                                   ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/index.rst
                                   ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/sdg/*
                                   ${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/sdg/sdglib/*
                                   ${DOXYGEN_INDEX_FILE}
                                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                                   COMMENT "Generating documentation with Sphinx, WD: ${CMAKE_CURRENT_BINARY_DIR}")
                add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})
            endif(PY_BREATHE)
        endif(SPHINX_FOUND AND NOT PYTHON_EXECUTABLE EQUAL "")
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
ENDIF()

IF(BUILD_PYTHON_INTERFACE)
    FIND_PACKAGE(PythonInterp 3 REQUIRED)
    FIND_PACKAGE(PythonLibs 3 REQUIRED)

    MESSAGE("Building python interface")
    FIND_PACKAGE(SWIG REQUIRED)
    INCLUDE(${SWIG_USE_FILE})
    set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/pysdg)
    configure_file(interfaces/python/__init__.py.in pysdg/__init__.py)
    configure_file(interfaces/SequenceDistanceGraph.i.in SequenceDistanceGraph.i)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTY SWIG_MODULE_NAME pysdg)

    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTIES CPLUSPLUS ON)
    if (NOT APPLE)
        set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTY COMPILE_DEFINITIONS "SWIGWORDSIZE64")
        set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTY GENERATED_COMPILE_DEFINITIONS "SWIGWORDSIZE64")
    endif()

    SET_PROPERTY(SOURCE ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTY CMAKE_SWIG_FLAGS "-doxygen")
    SWIG_ADD_LIBRARY(pysdg
                     LANGUAGE python
                     SOURCES ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i
                     )
    SWIG_LINK_LIBRARIES(pysdg ${PYTHON_LIBRARIES} sdg_static)
    set_target_properties(pysdg PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY pysdg)
    install(TARGETS
            pysdg
            RUNTIME DESTINATION pysdg
            ARCHIVE DESTINATION pysdg
            LIBRARY DESTINATION pysdg)
ENDIF()

IF(BUILD_R_INTERFACE)
    MESSAGE("Building R interface")
    FIND_PACKAGE(SWIG REQUIRED)
    SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
    FIND_PACKAGE(R REQUIRED)
    INCLUDE_DIRECTORIES(${R_INCLUDE_DIR})
    INCLUDE(${SWIG_USE_FILE})
    set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/rsdg)
    configure_file(interfaces/SequenceDistanceGraph.i.in SequenceDistanceGraph.i)

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTIES CPLUSPLUS ON)
    set_property(SOURCE ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i PROPERTY SWIG_MODULE_NAME rsdg)
    SWIG_ADD_LIBRARY(rsdg LANGUAGE r SOURCES ${CMAKE_BINARY_DIR}/SequenceDistanceGraph.i)
    SWIG_LINK_LIBRARIES(rsdg ${R_LIBRARIES} sdg_static)

    install(TARGETS
            rsdg
            RUNTIME DESTINATION rsdg
            ARCHIVE DESTINATION rsdg
            LIBRARY DESTINATION rsdg)
ENDIF()

install(TARGETS
        sdg-dbg
        sdg-workspace
        sdg-datastore
        sdg-kmercounts
        sdg-mapper
        sdg
        sdg_static
#        runBasicTests
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(DIRECTORY src/sdglib
        DESTINATION include)

IF (BUILD_JULIA_INTERFACE)
    message("Building Julia wrapper")
    #message("Using cmake prefix path: ${JlCxx_DIR}/lib/cmake/JlCxx")
    #message("Using jlcxx include: ${JlCxx_DIR}/include")
    #include_directories(${JlCxx_DIR}/include)
    find_package(JlCxx REQUIRED)
    get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
    get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
    message(STATUS "Founc JlCxx at ${JlCxx_location}")
        add_library(sdgjl SHARED interfaces/JuliaSDG.cpp)
        target_link_libraries(sdgjl sdg)
        target_link_libraries(sdgjl JlCxx::cxxwrap_julia)
    install(TARGETS sdgjl LIBRARY DESTINATION sdgjl ARCHIVE DESTINATION sdgjl RUNTIME DESTINATION sdgjl)
ENDIF()


