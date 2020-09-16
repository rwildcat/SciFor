# -----------------------------------------------------------------------------
# CMake for main src branch
# Ramon Solano <ramon.solano@gmail.com>
#
# apr/24/2020
#    initial version
# -----------------------------------------------------------------------------

# -- prerequisites
add_subdirectory (odepack)

# -- main target
add_library (scifor SciFor.f90 scifor_ode.f90)

# add odepack library (objects!)
target_link_libraries (scifor odepack)

# or include odepack library objects directly w/o library
# add_library(scifor SciFor.f90 scifor_ode.f90 $<TARGET_OBJECTS:odepack>)

# -- force to use preprocessing (even *.f90 instead *.F90)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")


# -- find_pkg for pc config file
find_package (PkgConfig)
if (PkgConfig_FOUND)
	# get path to pkg-config config files
	pkg_get_variable(PC_PATHS pkg-config pc_path)
	
	# make a list and get path 0
	string(REPLACE ":" ";" PC_PATH "${PC_PATHS}")
	list(GET PC_PATH 0 PC_PATH0)

	# init my CFLAGS and LDFLAGS
	set (MyCFLAGS "-I${CMAKE_INSTALL_PREFIX}/include")
	set (MyLDFLAGS "-L${CMAKE_INSTALL_PREFIX}/lib;-lscifor")
endif()


# -- FGSL
find_package (FGSL)
if (FGSL_FOUND)

	include_directories (${FGSL_INCLUDE_DIRS})
	target_link_libraries (scifor ${FGSL_LIBRARIES})
	add_compile_definitions(FGSL)
	
	# pk-config file
	if (PkgConfig_FOUND)
		pkg_search_module (FGSL_M fgsl)
		list (APPEND MyCFLAGS  ${FGSL_M_CFLAGS})
		list (APPEND MyLDFLAGS ${FGSL_M_LDFLAGS})
	endif()
	
else()
	message ("** WARNING: GNU Fortran Scientific Library (fgsl) not found")
endif()


# -- lapack
# macOS: macport install *pc files in /opt/local/lib/lapack/pkgconfig/
#	make sure to symlink or add to $PKG_CONFIG_PATH
find_package (LAPACK)
if (LAPACK_FOUND)

	message ("--  LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
	target_link_libraries (scifor ${LAPACK_LIBRARIES})
	add_compile_definitions (LAPACK)

	# pk-config file
	if (PkgConfig_FOUND)
		pkg_search_module (LAPACK_M lapack)
		if (LAPACK_M_MODULE_NAME)
			list (APPEND MyCFLAGS  ${LAPACK_M_CFLAGS})
			list (APPEND MyLDFLAGS ${LAPACK_M_LDFLAGS})
			message(" -->> MyLDFLAGS       =${MyLDFLAGS}")
		else()
			message ("-- WARNING: LAPACK pkg-config info NOT found")
		endif()
	endif()
else()
	message ("** WARNING LAPACK not found")

endif()

# -- FFTW3
find_package (FFTW3)
if (FFTW3_FOUND)

	include_directories (${FFTW3_INCLUDE_DIRS})
	target_link_libraries (scifor ${FFTW3_LIBRARIES})
	add_compile_definitions (FFTW)

	# pk-config file
	if (PkgConfig_FOUND)
		pkg_search_module (FFTW3_M fftw3)
		list (APPEND MyCFLAGS  ${FFTW3_M_CFLAGS})
		list (APPEND MyLDFLAGS ${FFTW3_M_LDFLAGS})
	endif()
	
else()
	message ("** WARNING: FFTW3 library not found")
endif()


# -- netCDF (from spacecenter.dk)
set (NETCDF_F90 "YES")
find_package (NetCDF)
if (NetCDF_FOUND)

	include_directories (${NETCDF_INCLUDES})
	target_link_libraries (scifor ${NETCDF_LIBRARIES})
	add_compile_definitions(NETCDF)

		# pk-config file
		if (PkgConfig_FOUND)
		pkg_search_module (netCDF_M netcdf-fortran)
		list (APPEND MyCFLAGS  ${netCDF_M_CFLAGS})
		list (APPEND MyLDFLAGS ${netCDF_M_LDFLAGS})
	endif()

else()
	message ("-- WARNING: netCDF NOT found")
endif()


# -- generate pk-config file
if (PkgConfig_FOUND)
	message ("-- Generating pkg-config file scifor.pc")

	# reformat from list to blank-cat string
	list (REMOVE_DUPLICATES MyCFLAGS)
	list (REMOVE_DUPLICATES MyLDFLAGS)
	
	list (JOIN MyCFLAGS " " MyCFLAGS)
	list (JOIN MyLDFLAGS " " MyLDFLAGS)

	configure_file (scifor.pc.in scifor.pc)
endif()


# --- install
install (TARGETS scifor DESTINATION lib)
install (FILES ${PROJECT_BINARY_DIR}/src/scifor.mod DESTINATION include)

if (PkgConfig_FOUND)
	install (FILES ${PROJECT_BINARY_DIR}/src/scifor.pc DESTINATION ${PC_PATH0})
endif()
