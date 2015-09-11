# adapted from https://raw.githubusercontent.com/lindahua/light-matrix/master/cmake_modules/FindMKL.cmake

# Copyright (c) 2012, Dahua Lin
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# a simple cmake script to locate Intel Math Kernel Library

# This script tries to find MKL in the environment variable MKLROOT and if that
# does not exist, it tries to find mkl.h

# Stage 1: find the root directory

set(MKLROOT_PATH $ENV{MKLROOT})

if (NOT MKLROOT_PATH)
  find_path (MKLROOT_PATH include/mkl.h)
endif ()

# Stage 2: find include path and libraries

if (MKLROOT_PATH)
  set(EXPECT_MKL_INCPATH "${MKLROOT_PATH}/include")

  if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib")
  endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")

  if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    if (CMAKE_SIZEOF_VOID_P MATCHES 8)
      set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib/intel64")
    else ()
      set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib/ia32")
    endif ()
  endif (CMAKE_SYSTEM_NAME MATCHES "Linux")

  # set include

  if (IS_DIRECTORY ${EXPECT_MKL_INCPATH})
    set(MKL_INCLUDE_DIRS ${EXPECT_MKL_INCPATH})
  endif ()

  if (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
  set(MKL_LIBRARY_DIR ${EXPECT_MKL_LIBPATH})
  endif ()

  # find specific library files
  if (CMAKE_SIZEOF_VOID_P MATCHES 8)
    find_library (LIB_MKL_INTEL NAMES mkl_intel_lp64 HINTS ${MKL_LIBRARY_DIR})
  else ()
    find_library (LIB_MKL_INTEL NAMES mkl_intel HINTS ${MKL_LIBRARY_DIR})
  endif ()
  find_library (LIB_MKL_CORE NAMES mkl_core HINTS ${MKL_LIBRARY_DIR})
  find_library (LIB_MKL_SEQUENTIAL NAMES mkl_sequential HINTS ${MKL_LIBRARY_DIR})

endif ()

set(MKL_LIBRARIES ${LIB_MKL_INTEL} ${LIB_MKL_CORE} ${LIB_MKL_SEQUENTIAL})

# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MKL DEFAULT_MSG
  MKL_LIBRARY_DIR
  MKL_LIBRARIES
  MKL_INCLUDE_DIRS)

mark_as_advanced(LIB_MKL_CORE LIB_MKL_INTEL LIB_MKL_SEQUENTIAL MKL_INCLUDE_DIRS)

if (CMAKE_SIZEOF_VOID_P MATCHES 8)
  set(MKL_STATIC_LIBRARIES "-Wl,-Bstatic,--start-group;mkl_intel_lp64;mkl_core;mkl_sequential;-Wl,--end-group,-Bdynamic")
else ()
  set(MKL_STATIC_LIBRARIES "-Wl,-Bstatic,--start-group;mkl_inte;mkl_core;mkl_sequential;-Wl,--end-group,-Bdynamic")
endif ()
