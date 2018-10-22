
#
# - This module locates an installed R distribution.
#
# Defines the following:
#  R_COMMAND           - Path to R command
#  R_HOME              - Path to 'R home', as reported by R
#  R_INCLUDE_DIR       - Path to R include directory
#  R_LIBRARY_BASE      - Path to R library
#  R_LIBRARY_BLAS      - Path to Rblas / blas library
#  R_LIBRARY_LAPACK    - Path to Rlapack / lapack library
#  R_LIBRARY_READLINE  - Path to readline library
#  R_LIBRARIES         - Array of: R_LIBRARY_BASE, R_LIBRARY_BLAS, R_LIBRARY_LAPACK, R_LIBRARY_BASE [, R_LIBRARY_READLINE]
#
# Variable search order:
#   1. Attempt to locate and set R_COMMAND
#     - If unsuccessful, generate error and prompt user to manually set R_COMMAND
#   2. Use R_COMMAND to set R_HOME
#   3. Locate other libraries in the priority:
#     1. Within a user-built instance of R at R_HOME
#     2. Within an installed instance of R
#     3. Within external system libraries
#

set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")
find_program(R_COMMAND R ../doc "R executable.")
set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

if(R_COMMAND)
    execute_process(WORKING_DIRECTORY ..
            COMMAND ${R_COMMAND} RHOME
            OUTPUT_VARIABLE R_ROOT_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    # the following command does nothing currently, but will be used when deprecated code is removed
    set(R_HOME ${R_ROOT_DIR} CACHE PATH "R home directory obtained from R RHOME")
    message("Found R: ${R_HOME}")

    find_path(R_INCLUDE_DIR R.h
            HINTS ${R_ROOT_DIR}
            PATHS /usr/local/lib /usr/local/lib64 /usr/share
            PATH_SUFFIXES include R/include
            )
    find_path(R_INCLUDE_DIR Rdefines.h
            HINTS ${R_ROOT_DIR}
            PATHS /usr/local/lib /usr/local/lib64 /usr/share
            PATH_SUFFIXES include R/include
            )

    MESSAGE("R include dir: ${R_INCLUDE_DIR}")

    find_library(R_LIBRARY_BASE R
            HINTS ${R_ROOT_DIR}/lib
            )

    find_library(R_LIBRARY_BLAS NAMES Rblas blas
            HINTS ${R_ROOT_DIR}/lib
            )

    find_library(R_LIBRARY_LAPACK NAMES Rlapack lapack
            HINTS ${R_ROOT_DIR}/lib
            )

    find_library(R_LIBRARY_READLINE readline
            )

else()
    message(SEND_ERROR "FindR.cmake requires the following variables to be set: R_COMMAND")
endif()

# Note: R_LIBRARY_BASE is added to R_LIBRARIES twice; this may be due to circular linking dependencies; needs further investigation
set(R_LIBRARIES ${R_LIBRARY_BASE} ${R_LIBRARY_BLAS} ${R_LIBRARY_LAPACK} ${R_LIBRARY_BASE})
message("R libraries path: ${R_LIBRARIES}")
if(R_LIBRARY_READLINE)
    set(R_LIBRARIES ${R_LIBRARIES} ${R_LIBRARY_READLINE})
endif()
