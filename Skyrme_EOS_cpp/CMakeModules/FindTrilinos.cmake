# adapted from https://raw.githubusercontent.com/dealii/dealii/master/cmake/modules/FindTRILINOS.cmake

## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the Trilinos library
#
# This module exports:
#
#   Trilinos_DIR
#   Trilinos_INCLUDE_DIRS
#   Trilinos_LIBRARIES
#

#
# If 'variable' is empty it will be set to 'value'
#
MACRO(SET_IF_EMPTY _variable)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${ARGN})
  ENDIF()
ENDMACRO()

SET(Trilinos_DIR "" CACHE PATH "An optional hint to a Trilinos installation")
SET_IF_EMPTY(Trilinos_DIR "$ENV{Trilinos_DIR}")

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(TRILINOS_CONFIG
  CONFIG QUIET
  NAMES Trilinos TRILINOS
  HINTS
    ${Trilinos_DIR}/lib/cmake/Trilinos
    ${Trilinos_DIR}
  PATH_SUFFIXES
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
    lib${LIB_SUFFIX}/cmake/Trilinos
  NO_SYSTEM_ENVIRONMENT_PATH
  )

#
# *Boy* Sanitize the include paths given by TrilinosConfig.cmake...
#
STRING(REGEX REPLACE
  "(lib64|lib)\\/cmake\\/Trilinos\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
  Trilinos_INCLUDE_DIRS "${Trilinos_INCLUDE_DIRS}"
  )

#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
SET(_libraries "")
FOREACH(_library ${Trilinos_LIBRARIES})
  LIST(APPEND _libraries Trilinos_LIBRARY_${_library})
  FIND_LIBRARY(Trilinos_LIBRARY_${_library}
    NAMES ${_library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )
ENDFOREACH()


find_package_handle_standard_args(Trilinos
  FOUND_VAR Trilinos_FOUND
  REQUIRED_VARS ${_libraries} Trilinos_TPL_LIBRARIES Trilinos_INCLUDE_DIRS
)
