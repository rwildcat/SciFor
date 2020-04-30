# `CMake` auxiliary tools

*R. Solano (<ramon.solano@gmail.com>)*   
*Last update: Apr/24/2020*

&nbsp;

**Contents:**

 `CMake` modules for extending default capabilities, such as for searching for additional packages, programs and libraries. 

For enabling CMake to use these modules, add the following code to the main `CMakeLists.txt` configuration file:

	set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

Once this extra modules are set in the `cmake` default modules path, these can be used as for example:

Module: 

	# find FGSL using my custom module `FindFGSL.cmake`
	find_package (FGSL)
	if (FGSL_FOUND)
		include_directories(${FGSL_INCLUDE_DIR})
		link_directories(${FGSL_LINK_DIRECTORIES})
		add_compile_definitions(FGSL)
	else(FGSL_FOUND)
		message(WARNING "** GNU Scientific Library (fgsl) NOT found **")
	endif(FGSL_FOUND)

		
`CMake` comes with various modules for finding various well-known libraries and packages. You can get a listing of which modules your version of CMake supports by typing 

	cmake --help-module-list

## Included custom modules

- **FindFGSL.cmake**: Finds the *GNU Fortran Scientific library*
	- References:  [CMake 3.17.1](https://cmake.org/cmake/help/v3.17/release/3.17.html#modules)
- **FindFFTW.cmake**: Finds the [*FFTW-3 library*](http://www.fftw.org)
	- Reference: [Hackage: The Haskell Package Repository](http://hackage.haskell.org/package/eigen-2.1.4/src/eigen3/cmake/FindFFTW.cmake)

.