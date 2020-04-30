# ----------------------------------------------------------
# SciFor - A modern Fortran interface to classic Fortran libraries
# (C) 2020 Ramon Solano <ramon.solano@gmail.com>
#
# Usage:
#	cd <proj_dir>; mkdir build; cd build
#	[set] cmake [opts] ..
#	make
#	mke install
#
#	cmake  [-DCMAKE_INSTALL_PREFIX=<path>] <main_dir>
#
# where some opts:
#	-DCMAKE_INSTALL_PREFIX = install dir (default: /usr/local)
#
# and some sets:
#	FC=<ccompiler>, e.g. FC=GFORTRAN
#
# For other available options:
# 	cmake --help
#
# Examples:
# 
# 1. Regular system-wide install
#	$ cd build
#	$ cmake ..
#	$ make
#	$ sudo make install
#
# 2. Personal installation  e.g. on ~/lib:
#	$ cd build
#	$ cmake -DCMAKE_INSTALL_PREFIX=~/lib ..
#	$ make
#	$ make install
#
# 3. To install in a different directory after the system was compiled:
#	 $ cd build
#   $ make install DESTDIR=<path>
#
# Optiona dependencies:
#	- FGSL / GSL
# 	- LAPACK
#  - FFTW
#  - NetCDF
# ----------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

enable_language(Fortran)
project(SciFor Fortran)

set (PACKAGE xbtm)
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (VERSION_PATCH 0)

# custom modules
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# main source
add_subdirectory (src)
