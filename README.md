SciFor
=====

A project to implement a modern Fortran API to specialized libs
---


> *Bringing high quality, extensively-proven classic math and engineering libraries to modern use* -- R. Solano

(C) 2020 - Ramon Solano <<ramon.solano@gmail.com>>

### --- *work in progress* ---


Fortran has been the prime computer programming language for scientific and engeneering purposes since the beginning of the computer era. A vast amount of highly-proven and curated mathematical algoritms has been developed using Fortran, and many of them are still in current usage.

As one of the key Fortran features has been backwards compatibility with older versions such as FORTRAN77, valuable algorithms written in old Fortran versions are kept in current usage. As an example, ODEPACK, from the Center for Applied Scientific Computing at the Lawrence Livermore National Laboratory. ODEPACK is a set of solvers for ODE Initial Value Problems. There are many more.

Unfortunately, original code can be challenging for new users due to its required preamble, especially for users accustomed to the flexible interfaces of modern languages, or it may be just a desirable feature to be able to use those powerful classic libraries using a modern, even an object oriented, interface.

Lets compare the syntax of the `LSODAR()` subroutine and as provided by `SciFor`s  `odesolv%solve()`:

* Original

	```fortran
	CALL DLSODAR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT, G, NG, JROOT)
	```
	
* SciFor

	```fortran	
	! example 1: t range 
	call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3] )
	
	! example 2: t range, root finding
	call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3], froots=fgex, nfroots=1 )
	
	! example 3: t points, rel/abs tolerances, jacobian, roots
	call ode%solve(fchem3, [1d0, 0d0, 0d0], tpts=[0d0, 4d-1, 4d0, 4d1, 4d2, 4d3, 4d4, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10], rtol=[1d-4], atol=[1.d-6, 1.d-10, 1.d-6], jac=fjchem3, froots=frchem3, nfroots=2 )
	```	
	
This project is aimed to promote the use of a formal lenaguaje such as Fortran among new users coming other modern languajes such as `Python`, etc., who may be demoralised for the steep, learning curve of the diverse classic Fortran libraries available. Or just to make life a bit easier to regular Fortran users.
	
Contributions are invited and welcome.
	
*Ramon Solano  
Colima, México.  
April/2020* 

### Libraries supported (basic functions so far)

* General      -- Pretty matrices printing
* NumFor       -- Numeric linear spaces, array flipping
* [LAPACK][]   -- Solve `Ax=b`
* [FGSL][]     -- Interpolation
* [FFTW][]     -- FFT and IFFT
* [NetCDF][]   -- Read data, attributes
* [ODEPACK][]  -- Solve sysems of ODEs (`LSODAR()`)
* [RANDLIB90][]  -- Random numbers

**NOTE**:

* The actual support is provided **according to the libraries effectively installed on the working computer**. For example, if your system only have  the [LAPACK][LAPACK] library, **only the LAPACK section will be enables**.

* There is a set of internal, own utilities such as `numfor_linspace()` which are independent of any other library. Those subprograms will be always available.

* The [ODEPACK][ODEPACK] library is included within, as packaged versions may be not available for your system. To date, it is the most recent version. If for any reason you want to use your local library, some modification to the `src/CMakeFile.txt` will be required. More on this in fiture versions.

Usage
-----

The easiest way is to use the provided `pkg-config` utility:

```sh
$ gfortran <your_program.f90> `pkg-config --cflags --libs scifor`
```

Otherwise, and depending upon your local libraries available:

```sh
$ gfortran <your_program.f90> -I/usr/local/include -I/opt/local/include/fgsl -I/opt/local/include -L/usr/local/lib -L/opt/local/lib -lscifor -framework Accelerate -lfgsl -lgsl -lgslcblas -lm -lfftw3 -lnetcdff
```


Installation
------------

The `SciFor` library is built using the `cmake` utility.

1. git clone
2. `$ cd SciFor`    
3. `$ cmake ..`    
4. `$ make`   
5. `$ sudo make install`


Examples
--------

Some usage examples are provided in the `./tests` directory. Excerpts from those examples are as follows:

### FFT (FFTW)

```fortran
program fft_scifor_test

	use scifor
	:
	real(8) :: y(N), yi(N)
	complex(8) :: yFFT(N)
	
	! make a composite signal
	y = ...

	! fft
	yFFT = fft_fft(y)  
	
	! your process here
	
	! ifft
	yi = fft_ifft(yFFT)
	
	:
```

### Interpolation (FGSL)

```fortran
program interp_scifor
	
	use scifor
	:
	real(8) :: xdata(ndat),  ydata(ndat)
	real(8) :: x(nplot), y(nplot)
	
	! field data
	xdata = [ 1.0, 2.1, 2.9, 3.8, 5.2, 6.4 ]
	ydata = [ 0.5, 3.4, 5.8, 4.1, 1.9, 0.5 ]
	
	! desired interp range
	call numfor_linspace(x, xdata(1), xdata(ndat))
	
	! linear, poly, cspline, akima, steffen, cspline_p, akima_p
	y = interp_interp(x, xdata, ydata, 'cspline')
	
	:
```

### NetCDF

```fortran
program netcdf_test

   use scifor
   :
   real, allocatable    :: tos(:,:,:)
   real(8), allocatable :: lat(:), lon(:)
   real                  :: fillval
 
   fname = 'tos_O1_2001-2002.nc'
   vname = 'tos' 
   xdim_name = 'lon'
   ydim_name = 'lat'
   fillname = '_FillValue'
   
   ! read lat, lon (1-d)
   call io_nc_getvar(fname, ydim_name, lat)
   call io_nc_getvar(fname, xdim_name, lon)

   ! get data (3-d)
   call io_nc_getvar(fname, vname, tos)
   
   ! get att fillvalue
   call io_nc_getatt(fname, vname, fillname, fillval)

   ! "mask" nodata, use min val
   where(var3d==fillval) var3d = minval(var3d)
   :
```

### Linear Algebra (LAPACK)

```fortran
program linalg1

   use scifor
	:   
   real(8) :: A(n,n)
   real(8) :: b(n), x(n)
   real(8) :: bn(n,nRHS), xn(n,nRHS)
   
   ! set data   
   A  = ...
   b  = ...
   bn = ...

   ! solve A*b = x
   x = linalg_solve(A, b)

   ! solve A*bn = xn
   xn = linalg_solve(A, bn)
   :
```

### ODE (ODEPACK -- LSODAR)

```fortran
program ode_scifor_test

   use scifor   
   type(Odesolv)  :: ode 
   
   ! t range
   call ode%solve ( fdy, [2d-3], trng=[0d0, 2d0/2d-3] )
   
   ! t range, root finding
   call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3], froots=fgex, nfroots=1 )
   
   ! ---

   ! example 2: van der Pol oscillator
   !  2nd deg, stiff/no stiff: 
   !  y'' = µ(1 - y^2)y' - y
   !  y(0), y'(0) = [2, 0]
   !  t = [0, 100]
   
   ! t range
   call ode%solve( fvdp, [2d0,0d0], trng=[0d0, 100d0] )
   call plt%plot(ode%t, ode%y(:,1))

   ! plus roots
   call ode%solve( fvdp, [2d0,0d0], trng=[0d0, 100d0], froots=frvdp, nfroots=1 )
   
   ! ---

   ! example 3: 3 chem substances ratios
   !  1st degree 3-eq system + t points + tolerances + jacobian + roots
   call ode%solve(                                             &
      fchem3, [1d0, 0d0, 0d0],                                 &
      tpts=[0d0, 4d-1, 4d0, 4d1, 4d2, 4d3, 4d4, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10],                                         &
      rtol=[1d-4], atol=[1.d-6, 1.d-10, 1.d-6],                &
      jac=fjchem3, froots=frchem3, nfroots=2 )
   :
```

### Random (RANLIB90)

```fortran
program test_random

   use scifor
   type(Scifor_random)  :: ran
   :
   real(4)              :: x, xv(n), xm(n,n)
   

   ! normal, beta, binomial, chisq, exp ...
   call ran%init('normal', [10.0, 1.0])
   call ran%rand(x)
   call ran%rand(xv)
   call ran%rand(xm)

   ! fixed seed n (n>1), e.g. 2
   call ran%init('normal', [10.0, 1.0], 2)
   call ran%rand(...)

   ! reset seed, keep params
   call ran%init(seed=2)
   :
```

[LAPACK]: http://www.netlib.org/lapack/
[FGSL]: https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library
[FFTW]: http://www.fftw.org
[NetCDF]: https://www.unidata.ucar.edu/software/netcdf/
[ODEPACK]: https://computing.llnl.gov/casc/odepack/
[RANDLIB90]: https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/27
[ford]: https://github.com/cmacmackin/ford
[Ford]: https://github.com/Fortran-FOSS-Programmers/ford/wiki
