SciFor
=====

A project to implement a modern Fortran API to specialized classic libs
---


> *Bringing high quality, extensively-proven classic math and engineering libraries to modern use* -- R. Solano

(C) 2020 - Ramon Solano <<ramon.solano@gmail.com>>

Fortran has been the prime computer programming language since the beginning of the scientific computing era. A vast amount of highly-proven and curated mathematical algoritms has been developed since the early versions, and many of them are still in current usage.

As one of the key Fortran features has been backwards compatibility with older versions such as FORTRAN77, valuable algorithms written in old Fortran versions are kept in current usage. As an example, ODEPACK, from the Center for Applied Scientific Computing at the Lawrence Livermore National Laboratory. ODEPACK is a set of solvers for ODE Initial Value Problems.

Lets compare the syntax of the `LSODAR()` subroutine:

* Original

	```fortran
	DLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT,
     2            G, NG, JROOT)
	```
	
* SciFor

	```fortran
	type(Odesolv) :: ode 
	
	! example 1: 
	!	1st deg, stiff + roots: flame radius, y' = y^2 - y^3
	call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3] )
	
	! example 2:
	!	2nd deg, stiff/no stiff: van der Pool oscillator, y'' = m(1 - y^2)y' - y
	call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3], froots=fgex, nfroots=1 )
	
	! example 3
	!	1st deg, 3-eq system + roots + jacobian: 3 substances ratios
	call ode%solve(                &
       fchem3, [1d0, 0d0, 0d0],    &
       tpts=[0d0, 4d-1, 4d0, 4d1, 4d2, 4d3, 4d4, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10],  &
       rtol=[1d-4], atol=[1.d-6, 1.d-10, 1.d-6],   &
       jac=fjchem3, froots=frchem3, nfroots=2 )
	```	
	
This project is aimed to promote the use of a formal lenaguaje such as Fortran among new users coming other modern languajes such as `Python`, etc., who may be demoralised for the steep, sometimes heteregeneous, learning curve of the diverse classic Fortran libraries available.
	
Contributions are invited and welcome.
	
*Ramon Solano  
Colima, México.  
April/2020* 

### Libraries supported (basic functions so far)

* [LAPACK][LAPACK]   -- `linalg_solve()`
* [FGSL][FGSL]     -- `interp()`
* [FFTW][FFTW]     -- `fft()`, `ifft()`
* [NetCDF][NetCDF]   -- `get_var()`, `get_att()`, `getvardims()`
* [ODEPACK][ODEPACK]  --  `LSODAR()` 

**NOTE**:

* The actual support is provided **according to the libraries effectively installed on the working computer**. For example, if your system only have  the [LAPACK][LAPACK] library, **only the LAPACK section will be enables**.
* There is a set of internal, own utilities such as `numfor_linspace()` which are independent of any other library. Those subprograms will be always available.
* The [ODEPACK][ODEPACK] library is included within, as packaged versions may be not available for your system. To date, it is the most recent version. If for any reason you want to use your local library, some modification to the `src/CMakeFile.txt` will be required. More on this in fiture versions.

Installation
------------

The `SciFor` library is built using the `cmake` utility.

1.Git clone   
2. `$ cd SciFor`    
3. `$ cmake ..`    
4. `$ make`   
5. `$ sudo make install`




[LAPACK]: http://www.netlib.org/lapack/
[FGSL]: https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library
[FFTW]: http://www.fftw.org
[NetCDF]: https://www.unidata.ucar.edu/software/netcdf/
[ODEPACK]: https://computing.llnl.gov/casc/odepack/
[Ford]: https://github.com/Fortran-FOSS-Programmers/ford/wiki