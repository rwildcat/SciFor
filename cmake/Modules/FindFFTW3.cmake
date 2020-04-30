# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindFFTW3
--------

Find the FFTW3 includes and libraries.

FFTW3 is a C subroutine library for computing the discrete Fourier transform (DFT) in one or more dimensions, of arbitrary input size, and of both real and complex data (as well as of even/odd data, i.e. the discrete cosine/sine transforms or DCT/DST). 

Imported Targets
^^^^^^^^^^^^^^^^

If FFTW3 is found, this module defines the following :prop_tgt:`IMPORTED`
targets::

 FFTW3::fftw3      - The main FFTW3 library.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project::

 FFTW3_FOUND          - True if FFTW3 found on the local system
 FFTW3_INCLUDE_DIRS   - Location of FFTW3 header files.
 FFTW3_LIBRARIES      - The FFTW3 libraries.
 FFTW3_VERSION        - The version of the discovered FFTW3 install.

Hints
^^^^^

Set ``FFTW3_ROOT_DIR`` to a directory that contains a FFTW3 installation.

This script expects to find libraries at ``$FFTW3_ROOT_DIR/lib`` and the FFTW3
headers at ``$FFTW3_ROOT_DIR/include``.  The library directory may
optionally provide Release and Debug folders. If available, the library
named ``fftw3d``, is recognized as debug library.

Cache Variables
^^^^^^^^^^^^^^^

This module may set the following variables depending on platform and type
of FFTW3 installation discovered.  These variables may optionally be set to
help this module find the correct files::

 FFTW3_LIBRARY             - Location of the FFTW3 library.
 FFTW3_LIBRARY_DEBUG       - Location of the debug FFTW3 library (if any).

#]=======================================================================]

include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``FFTW3_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{FFTW3_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{FFTW3_ROOT_DIR}" FFTW3_ROOT_DIR )
  set( FFTW3_ROOT_DIR "${FFTW3_ROOT_DIR}" CACHE PATH "Prefix for FFTW3 installation." )
endif()

if( NOT EXISTS "${FFTW3_ROOT_DIR}" )
  set( FFTW3_USE_PKGCONFIG ON )
endif()

#=============================================================================
# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``FFTW3_INCLUDEDIR`` and ``FFTW3_LIBDIR`` used below.
if( FFTW3_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( FFTW3 QUIET fftw3 )

  if( EXISTS "${FFTW3_INCLUDEDIR}" )
    get_filename_component( FFTW3_ROOT_DIR "${FFTW3_INCLUDEDIR}" DIRECTORY CACHE)    
  endif()
endif()

#=============================================================================
# Set FFTW3_INCLUDE_DIRS and FFTW3_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $FFTW3_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ${FFTW3_ROOT_DIR}/include ${FFTW3_INCLUDEDIR}
)

find_library( FFTW3_LIBRARY
  NAMES fftw3
  HINTS ${FFTW3_ROOT_DIR}/lib ${FFTW3_LIBDIR}
  PATH_SUFFIXES Release Debug
)

# Do we also have debug versions?
find_library( FFTW3_LIBRARY_DEBUG
  NAMES fftw3-dbg
  HINTS ${FFTW3_ROOT_DIR}/lib ${FFTW3_LIBDIR}
  PATH_SUFFIXES Debug
)

set( FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR} )
set( FFTW3_LIBRARIES ${FFTW3_LIBRARY})


#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( FFTW3
  FOUND_VAR
    FFTW3_FOUND
  REQUIRED_VARS
    FFTW3_INCLUDE_DIR
    FFTW3_LIBRARY
  VERSION_VAR
    FFTW3_VERSION
    )

mark_as_advanced( FFTW3_ROOT_DIR  FFTW3_LIBRARY FFTW3_INCLUDE_DIR
   FFTW3_LIBRARY_DEBUG FFTW3_USE_PKGCONFIG FFTW3_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" FFTW3_LIBRARY_DLL       "${FFTW3_LIBRARY}" )
  string( REPLACE ".lib" ".dll" FFTW3_LIBRARY_DEBUG_DLL "${FFTW3_LIBRARY_DEBUG}" )
endif()

if( FFTW3_FOUND AND NOT TARGET FFTW3::fftw3 )
  if( EXISTS "${FFTW3_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( FFTW3::fftw3      SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( FFTW3::fftw3 PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${FFTW3_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FFTW3_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "f90" )

    # If we have both Debug and Release libraries
    if( EXISTS "${FFTW3_LIBRARY_DEBUG_DLL}")
      set_property( TARGET FFTW3::fftw3 APPEND PROPERTY IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( FFTW3::fftw3 PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${FFTW3_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${FFTW3_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create
    # the imported library targets.
    add_library( FFTW3::fftw3      UNKNOWN IMPORTED )
    set_target_properties( FFTW3::fftw3 PROPERTIES
      IMPORTED_LOCATION                 "${FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FFTW3_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "f90" )
  endif()
endif()
