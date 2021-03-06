# -----------------------------------------------------------------------------
# CMake for main src branch
# Ramon Solano <ramon.solano@gmail.com>
#
# apr/24/2020
#		initial version
# sep/16/20
#		use `-framework Accelerate` If -DBLA_VENDOR undefined and oan Mac OS
# -----------------------------------------------------------------------------

# -- prerequisites
add_subdirectory (odepack)
add_subdirectory (randlib)

# -- main target
# add local libraries (as included objects!)
add_library (scifor SciFor.f90 
		scifor_ode.f90 scifor_random.f90
		$<TARGET_OBJECTS:odepack>
		$<TARGET_OBJECTS:randlib>  
		)

# included directories
target_include_directories(scifor PUBLIC ${CMAKE_BINARY_DIR}/src/randlib)

# -- force to use preprocessing for DEFINED statement
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# -- find_pkg for pc config file
find_package (PkgConfig)
if (PkgConfig_FOUND)
	# get path to pkg-config config files
	pkg_get_variable(PC_PATHS pkg-config pc_path)
	
	# make a list and get path 0
	string(REPLACE ":" ";" PC_PATH "${PC_PATHS}")
	list(GET PC_PATH 0 PC_PATH0)

	# init my CFLAGS and LDFLAGS for scifor .pc file
	set (MyCFLAGS "-I${CMAKE_INSTALL_PREFIX}/include")
	set (MyLDFLAGS "-L${CMAKE_INSTALL_PREFIX}/lib;-lscifor")
endif()

# -- lapack
# macOS: macport install *pc files in /opt/local/lib/lapack/pkgconfig/
#	make sure to symlink or add to $PKG_CONFIG_PATH

# warn about the `Accelerated` framework
if (APPLE AND NOT BLA_VENDOR)
	message("** Mac OS detected. Using Apples 'Accelerate' framework")
	message("** To override, use: `cmake .. -DBLA_VENDOR=<XXX>`; XXX=All, Intel, etc.")
	set (BLA_VENDOR "Apple")
endif()

find_package (LAPACK)
if (LAPACK_FOUND)
	message(>> "LAPACK found: LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
	target_link_libraries (scifor ${LAPACK_LIBRARIES})
	add_compile_definitions (LAPACK)

	# pk-config info file
	if (BLA_VENDOR STREQUAL "Apple")
			list (APPEND MyLDFLAGS "-framework Accelerate")
	else()
		if (PkgConfig_FOUND)
				pkg_search_module (LAPACK_M lapack)
				if (LAPACK_M_MODULE_NAME)
					list (APPEND MyCFLAGS  ${LAPACK_M_CFLAGS})
					list (APPEND MyLDFLAGS ${LAPACK_M_LDFLAGS})					
			else()
				message ("** WARNING: LAPACK pkg-config info NOT found")
			endif()
		endif()
	endif()

else()
	message ("** WARNING LAPACK not found")
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
