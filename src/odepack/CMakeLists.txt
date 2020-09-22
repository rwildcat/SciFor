# -----------------------------------------------------------------------------
# CMake for main src branch
# Ramon Solano <ramon.solano@gmail.com>
#
# apr/24/2020
#    initial version
# -----------------------------------------------------------------------------


# -- main target
add_library (odepack OBJECT opkdmain.f opkda1.f opkda2.f)

# -- set fortran dialect to legacy
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy -w")
