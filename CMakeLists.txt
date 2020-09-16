# --------------------------------------------------------------------
# SciFor - A modern Fortran interface to classic Fortran libraries
# (C) 2020 Ramon Solano <ramon.solano@gmail.com>
#
# Optional dependencies:
#	- FGSL / GSL
# 	- LAPACK
#  - FFTW
#  - NetCDF
#  - MINPACK
# --------------------------------------------------------------------

cmake_minimum_required (VERSION 3.0)

project (scifor 
			VERSION 0.2.0
			DESCRIPTION "SciFor - A modern Fortran interface to classic Fortran libraries"
			LANGUAGES Fortran)

# custom modules
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# main source
add_subdirectory (src)
