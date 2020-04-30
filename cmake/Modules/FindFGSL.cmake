# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindFGSL
--------

Find the native GNU Fortran Scientific Library (FGSL) includes and libraries.

The GNU Fortran Scientific Library (FGSL) is a numerical library for Fortran
programmers. It is free software under the GNU General Public
License.

Imported Targets
^^^^^^^^^^^^^^^^

If FGSL is found, this module defines the following :prop_tgt:`IMPORTED`
targets::

 FGSL::fgsl      - The main FGSL library.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project::

 FGSL_FOUND          - True if FGSL found on the local system
 FGSL_INCLUDE_DIRS   - Location of FGSL header files.
 FGSL_LIBRARIES      - The FGSL libraries.
 FGSL_VERSION        - The version of the discovered FGSL install.

Hints
^^^^^

Set ``FGSL_ROOT_DIR`` to a directory that contains a FGSL installation.

This script expects to find libraries at ``$FGSL_ROOT_DIR/lib`` and the FGSL
headers at ``$FGSL_ROOT_DIR/include/fgsl``.  The library directory may
optionally provide Release and Debug folders. If available, the libraries
named ``fgsld``, ``fgslblasd`` or ``cblasd`` are recognized as debug libraries.
For Unix-like systems, this script will use ``$FGSL_ROOT_DIR/bin/fgsl-config``
(if found) to aid in the discovery of FGSL.

Cache Variables
^^^^^^^^^^^^^^^

This module may set the following variables depending on platform and type
of FGSL installation discovered.  These variables may optionally be set to
help this module find the correct files::

 FGSL_CONFIG_EXECUTABLE   - Location of the ``fgsl-config`` script (if any).
 FGSL_LIBRARY             - Location of the FGSL library.
 FGSL_LIBRARY_DEBUG       - Location of the debug FGSL library (if any).

#]=======================================================================]

include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``FGSL_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{FGSL_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{FGSL_ROOT_DIR}" FGSL_ROOT_DIR )
  set( FGSL_ROOT_DIR "${FGSL_ROOT_DIR}" CACHE PATH "Prefix for FGSL installation." )
endif()

if( NOT EXISTS "${FGSL_ROOT_DIR}" )
  set( FGSL_USE_PKGCONFIG ON )
endif()

#=============================================================================
# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``FGSL_INCLUDEDIR`` and ``FGSL_LIBDIR`` used below.
if( FGSL_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( FGSL QUIET fgsl )

  if( EXISTS "${FGSL_INCLUDEDIR}" )
    get_filename_component( FGSL_ROOT_DIR "${FGSL_INCLUDEDIR}" DIRECTORY CACHE)    
  endif()
endif()

#=============================================================================
# Set FGSL_INCLUDE_DIRS and FGSL_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $FGSL_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( FGSL_INCLUDE_DIR
  NAMES fgsl.mod
  HINTS ${FGSL_ROOT_DIR}/include ${FGSL_INCLUDEDIR}
)

find_library( FGSL_LIBRARY
  NAMES fgsl
  HINTS ${FGSL_ROOT_DIR}/lib ${FGSL_LIBDIR}
  PATH_SUFFIXES Release Debug
)

# Do we also have debug versions?
find_library( FGSL_LIBRARY_DEBUG
  NAMES fgsld fgsl
  HINTS ${FGSL_ROOT_DIR}/lib ${FGSL_LIBDIR}
  PATH_SUFFIXES Debug
)

set( FGSL_INCLUDE_DIRS ${FGSL_INCLUDE_DIR} )
set( FGSL_LIBRARIES ${FGSL_LIBRARY})


#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set FGSL_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( FGSL
  FOUND_VAR
    FGSL_FOUND
  REQUIRED_VARS
    FGSL_INCLUDE_DIR
    FGSL_LIBRARY
  VERSION_VAR
    FGSL_VERSION
    )

mark_as_advanced( FGSL_ROOT_DIR  FGSL_LIBRARY FGSL_INCLUDE_DIR
   FGSL_LIBRARY_DEBUG FGSL_USE_PKGCONFIG FGSL_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" FGSL_LIBRARY_DLL       "${FGSL_LIBRARY}" )
  string( REPLACE ".lib" ".dll" FGSL_LIBRARY_DEBUG_DLL "${FGSL_LIBRARY_DEBUG}" )
endif()

if( FGSL_FOUND AND NOT TARGET FGSL::fgsl )
  if( EXISTS "${FGSL_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( FGSL::fgsl      SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( FGSL::fgsl PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${FGSL_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${FGSL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FGSL_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "f90" )

    # If we have both Debug and Release libraries
    if( EXISTS "${FGSL_LIBRARY_DEBUG_DLL}")
      set_property( TARGET FGSL::fgsl APPEND PROPERTY IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( FGSL::fgsl PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${FGSL_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${FGSL_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create
    # the imported library targets.
    add_library( FGSL::fgsl      UNKNOWN IMPORTED )
    set_target_properties( FGSL::fgsl PROPERTIES
      IMPORTED_LOCATION                 "${FGSL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FGSL_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "f90" )
  endif()
endif()
