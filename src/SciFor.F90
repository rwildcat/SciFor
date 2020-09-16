!># SciFor: A modern Fortran API to specialized libs
!
! Author: Ramon Solano
!			 ramon.solano@gmail.com
!			 Colima, Mexico
!
!## Contents
!
! 1. Linear algebra
!	Computes the solution to system of linear equations A * X = B for GE matrices
!	* SciFor	: linag_solve()
!	* Library: `LAPACK` (<http://www.netlib.org/lapack/>)
!	* Source : `dgesv()`
!
! 2. Interpolation
!	Interpolates a 1-D array \(y_i\) corresponding to a set of \(x_i\) values, 
! based on known \((x_j,y_j)\) point values. Includes Cubic, Akima, and Steffen
! splines algorithms.
!	* SciFor	 : interp_interp()
!  * Library : 
!		* [`FGSL`](https://doku.lrz.de/display/PUBLIC/FGSL+-+A+Fortran+interface+to+the+GNU+Scientific+Library)
!		* [`GSL`](https://www.gnu.org/software/gsl/doc/html/index.html)
!	* Source	 : [`interp()`](https://www.gnu.org/software/gsl/doc/html/interp.html)
!
! 3. FFT 
!	Computes the discrete Fourier transform in 1D arrays, of arbitrary input
! size, and of both real and complex data (as well as of even/odd data, i.e.
! the discrete cosine/sine transforms or DCT/DST).
!	* SciFor	 : `fft_fft()`, `fft_ifft()`
!	* Library : [`FFTW`](http://www.fftw.org)
!	* Source  : `dfftw_execute_dft_r2c()`
!
! 4. NetCDF
! 	Reads data and attributes from NetCDF files.
!	* SciFor  : `io_netcdf_getvar()`, `io_netcdf_getatt()`
!	* Library :  [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf
!	* Source  : `nf90_get_var()`, `nf90_get_att()`
!
! 5. ODEPACK
!	 Solves systems dy/dt = f with a dense or banded Jacobian when the problem is
! stiff. It automatically selects between nonstiff (Adams) and stiff (BDF) 
! methods. Additionally, a rootfinding capability is added.
!	* SciFor  : ode%solve()
!	* Library : ODEPACK
!	* Source  : (), which automatically switches between nonstiff and stiff
! solvers. Includes the solving for solutions of a related set of equations
! (roots).
!	* Source  : [`DLSODAR`](https://www.netlib.org/odepack/opkd-sum)
! 
!
!## To compile:
!
!- In-line: 
!
!```sh
!  $ gfortran SciFor.f90 xxx.f90 `pkg-config --cflags --libs fgsl fftw3 netcdf-fortran` -logpf -llapack
!```
!
!- As independent object:
!
!```sh
!  $ gfortran -c SciFor.f90  `pkg-config --cflags --libs fgsl fftw3 netcdf-fortran` -logpf -llapack
!```
!
!- As external library
!
!```sh
!  $ cd build
!  $ cmake ..
!  $ make
!  $ sudo make install
!```
!

module SciFor
   
   use scifor_ode
   implicit none

   ! default
   private

   ! except
   public   :: Odesolv
   public   :: io_print_mat, 	numfor_linspace
    				
#ifdef LAPACK    				 
	public	::  linalg_solve
#endif

#ifdef FFTW
	public	::  fft_fft, fft_ifft
#endif

#ifdef FGSL
	public 	:: interp_interp
#endif

#ifdef NETCDF
	public	:: io_nc_getvar, io_nc_getvardims, io_nc_getatt
#endif
					

   
   ! ---
   ! high-level interfaces 
   ! ---   

   !> Print a matrix with format, a row per matrix row.
   !
   !### Usage
   !
   !```fortran
   !call io_print_mat(A [, fmt])
   !```
   !
   ! Where:
   !
   ! `A`    : Matrix to print (`real(4)`, `real(8)`, `integer`)     
   ! `fmt`  : Print format, optional (defaults: 'ES12.4' for reals and 'I10' for integer)    
   !
   !### Example
   !
   !```fortran
   ! real(8) :: a(3,4)
   ! call io_print_mat(a)
   ! call io_print_mat(a, 'F10.4')
   !```
   interface io_print_mat
      module procedure print_mat_r8
      module procedure print_mat_r4
      module procedure print_mat_i
   end interface io_print_mat


   !> Flip an array
   !
   !### Usage
   !
   !```fortran
   !call numfor_flip(a)
   !```
   !
   ! Where:
   !
   ! `a` : 1-D array to flip (`real(4)`, `integer(4)`, `character`)     
   !
   !### Example
   !
   !```fortran
   !real(4) :: a(10)
   !character(40) :: c(5)
   !
   !call numfor_flip(a)
   !call numfor_flip(c)
   !```
   interface numfor_flip
         module procedure flip_1d_r4
         module procedure flip_1d_i4
         module procedure flip_1d_ch
   end interface
   
   !> Make a linespace between `[a,b]` (inclusive)
   !
   !### Usage
   !
   !```fortran
   !call numfor_linspace(x, a, b)
   !```
   !
   ! Where:
   !
   ! `x` : 1-D array where to store the lineal space (`real(8)`)     
   ! `a`, `b` : Interval limits, inclusive (`real(8)`, `real(4)`, `integer`)
   !
   !### Example
   !
   !```fortran
   ! real(8) :: x(100)
   ! call numfor_linspace(x, -6.28, 6.28)
   ! call numfor_linspace(x, 0, 10)
   !```
   interface numfor_linspace
      module procedure numfor_linspace8
      module procedure numfor_linspace4
      module procedure numfor_linspacei
   end interface

#ifdef FGSL
   !> Interpolation using [FGSL]() (Fortran [GNU Scientific Library]())
   !
   ! Interpolates a 1-D array \(y_i\) corresponding to a set of \(x_i\) values, based on 
   ! known \((x_j,y_j)\) point values.
   !
   ! **Interpolation methods:**
   !
   ! - `linear`   : Linear interpolation
   ! - `poly`     : Polynomial interpolation. This method should only be used for interpolating
   !    small numbers of points because polynomial interpolation introduces large oscillations, 
   !    even for well-behaved datasets.
   ! - `cspline`  : Cubic spline with natural boundary conditions. The resulting curve is piecewise 
   !    cubic on each interval, with matching first and second derivatives at the supplied data-points. 
   !    The second derivative is chosen to be zero at the first point and last point.
   ! - `akima`    : Non-rounded Akima spline with natural boundary conditions. This method uses the 
   !    non-rounded corner algorithm of Wodicka.
   ! - `steffen`   : Steffen’s method guarantees the monotonicity of the interpolating function between
   !    the given data points. Therefore, minima and maxima can only occur exactly at the data points, 
   !    and there can never be spurious oscillations between data points. The interpolated function is 
   !    piecewise cubic in each interval.
   ! - `cspline_p` : Cubic spline with *periodic* boundary conditions. The resulting curve is piecewise 
   !    cubic on each interval, with matching first and second derivatives at the supplied data-points. 
   !    The derivatives at the first and last points are also matched. Note that the last point in the 
   !    data must have the same y-value as the first point, otherwise the resulting periodic interpolation 
   !    will have a discontinuity at the boundary.
   ! - `akima_p`   : Non-rounded Akima spline with *periodic* boundary conditions. This method uses the 
   !    non-rounded corner algorithm of Wodicka.
   !
   ! The interpolation is piecewise smooth, and its behavior at the end-points is determined 
   ! by the type of interpolation used.
   !
   ! These interpolation methods produce curves that pass through each datapoint. 
   ! To interpolate noisy data with a smoothing curve see [[Splines]].
   !
   !Example
   !-------
   !
   !```fortran
	! integer, parameter :: ndat=6, nplot=100
   ! real(8) :: xdata(ndat),  ydata(ndat)
   ! real(8) :: xplot(nplot), yplot(nplot)
   !
	! ! field data
	! xdata = [ 1.0, 2.1, 2.9, 3.8, 5.2, 6.4 ]
	! ydata = [ 0.5, 3.4, 5.8, 4.1, 1.9, 0.5 ]
	!
	! !  interp range
	! call numfor_linspace(xplot, xdata(1), xdata(ndat))
	!
	! ! interpolate: linear, poly, cspline, akima, steffen, cspline_p, akima_p
   ! yplot = interp_interp1d(xplot, xdata, ydata, 'poly')
   !
   !```
   interface interp_interp
         module procedure interp_interp1d
   end interface
#endif

#ifdef LAPACK   
   !> Solve \(x\) in \(Ax = b\) (LAPACK)
   interface linalg_solve
      module procedure linalg_solve_b1    
      module procedure linalg_solve_bn
   end interface
#endif

#ifdef FFTW   
   !> fft ([[FFTW-3]])
   interface fft_fft
      module procedure fft_1d_r2c
   end interface
#endif

#ifdef FFTW   
   ! ifft (FFTW-3)
 	interface fft_ifft
      module procedure ifft_1d_c2r
   end interface
#endif
   
!  netcdf getvardims
   ! interface io_nc_getvardims
   !    module procedure io_nc_getvardims
   ! end interface

#ifdef NETCDF
   !> NetCDF get var 
   interface io_nc_getvar
      module procedure io_nc_getvar_1d_r8
      module procedure io_nc_getvar_2d_r4
      module procedure io_nc_getvar_3d_r4
   end interface
#endif

#ifdef NETCDF
   !> NetCDF get attribute
   interface io_nc_getatt
      module procedure io_nc_getatt_1d_r4
      module procedure io_nc_getatt_r4
   end interface
#endif

contains

   ! --- 
   ! low-level detailed subprograms 
   ! ---

   !> Flips a 1-d, r4 array
   subroutine flip_1d_r4(a)

      ! parameters
      real(4), intent(inout) :: a(:)   !! Array to flip
      
      ! local vars
      real(4), allocatable :: tmp(:)

      ! get mem
      allocate(tmp(size(a)))
      
      ! flip
      tmp = a
      a = tmp(size(a):1:-1)

      ! cleanup
      deallocate(tmp)

   end subroutine flip_1d_r4


   !> Flips a 1-d, i4 array
   subroutine flip_1d_i4(a)
      integer(4), intent(inout) :: a(:)
      ! local vars
      integer(4), allocatable :: tmp(:)

      ! get mem
      allocate(tmp(size(a)))
      
      ! flip
      tmp = a
      a = tmp(size(a):1:-1)

      ! cleanup
      deallocate(tmp)

   end subroutine flip_1d_i4


   !> Flips a 1-d, character array
   subroutine flip_1d_ch(a)
      character(*), intent(inout) :: a(:)
      ! local vars
      character(999) :: tmp
      integer :: n, i

      n = size(a)
      do i=1, n/2
         tmp = a(i)
         a(i) = a(n-i+1)
         a(n-i+1) = trim(tmp)
      end do

   end subroutine flip_1d_ch
   

   !> Makes a linespace between `[a,b]`. The number of steps is accordining `size(x)`.
   subroutine numfor_linspace8(x,a,b)
      real(8), intent(out) :: x(:)  !! 1-D array to save computed linespace into
      real(8), intent(in) :: a, b   !! Linespace limits (real(8), inclusive)
      
      ! my vars
      real(8) :: dx 
      integer :: i
      
      ! process
      dx = (b-a)/(size(x)-1)
      x = [(a+i*dx, i=0, size(x))]
      
   end subroutine numfor_linspace8
   
   ! ---
   
   !> Makes a linespace between `[a,b]`. The number of steps is accordining `size(x)`.
   subroutine numfor_linspace4(x,a,b)
      real(8), intent(out) :: x(:)  !! 1-D array to save computed linespace into
      real(4), intent(in) :: a, b   !! Linespace limits (real(4), inclusive)
      
      call numfor_linspace8(x, real(a,8), real(b,8))     
      
   end subroutine numfor_linspace4
     
   ! ---

   !> Makes a linespace between `[a,b]`. The number of steps is accordining `size(x)`.
   subroutine numfor_linspacei(x,a,b)  
      real(8), intent(out) :: x(:)    !! 1-D array to save computed linespace into
      integer, intent(in) :: a, b     !! Linespace limits (integer, inclusive)
      
      call numfor_linspace8(x, real(a,8), real(b,8))     
      
   end subroutine numfor_linspacei
   
   ! ---
    
   !> Print the real(8) matrix `A` with format `fmt_`, a row per matrix row
   subroutine print_mat_r8(A, fmt_)
      
      real(8), intent(in) :: A(:,:)                !! matrix to print
      character(*),optional, intent(in) :: fmt_    !! print format (default: 'ES12.4').
      
      ! my vars
      character(80) :: nc, fmt, fmt0
      integer :: i

      fmt = 'ES12.4'
      
      if (present(fmt_)) then
         fmt = fmt_
      end if

      write(nc,*) size(A,2)
      nc = adjustl(nc)
   
      write (fmt0, *)  '(' // trim(nc) // '(' // trim(fmt) // '))'
   
      do i=1, size(A,1)
         print fmt0, A(i,:)
      end do
      
   end subroutine print_mat_r8
    

   !> Print the real(4) matrix `A` with format `fmt_`, a row per matrix row
   subroutine print_mat_r4(A, fmt_)
      real(4), intent(in) :: A(:,:)                !! matrix to print
      character(*),optional, intent(in) :: fmt_    !! print format (default: 'ES12.4')
      
      if (present(fmt_)) then
         call print_mat_r8(real(A,8), fmt_)
      else
         call print_mat_r8(real(A,8))
      end if
      
   end subroutine print_mat_r4


   !> Print the integer matrix `A` with format `fmt_`, a row per matrix row
   subroutine print_mat_i(A, fmt_)
      integer, intent(in) :: A(:,:)                !! matrix to print
      character(*),optional, intent(in) :: fmt_    !! print format (default: 'I10')
      
      ! my vars
      character(80) :: nc, fmt, fmt0
      integer :: i

      fmt = 'I10'
      
      if (present(fmt_)) then
         fmt = fmt_
      end if

      write(nc,*) size(A,2)
      nc = adjustl(nc)
   
      write (fmt0, *)  '(' // trim(nc) // '(' // trim(fmt) // '))'
   
      do i=1, size(A,1)
         print fmt0, A(i,:)
      end do
      
   end subroutine print_mat_i

   ! ---

#ifdef LAPACK   
   ! Computes the solution to a real system of linear equations A * X = B,
   ! where :
   !     A = NxN coefficients matrix 
   !     B = RHS matrix (N x NRHS)
   !     X = Solution matrix (N x NRHS)
   !
   ! The LU decomposition with partial pivoting and row interchanges is
   ! used to factor A as A = P * L * U,
   ! where P is a permutation matrix, L is unit lower triangular, and U is
   ! upper triangular.  The factored form of A is then used to solve the
   ! system of equations A * X = B.
   !
   ! dSave_ = Save A and b (do not overwrite ) : default = True
   ! info_  = Return lapack exit code
   function linalg_solve_bn(A, b, dSave_, info_) result (x)
      real(8), intent (in out) :: a(:,:)
      real(8), intent (in out) :: b(:,:)
      logical, optional, intent(in)  :: dSave_
      integer, optional, intent(out) :: info_
      real(8) :: x(size(b,1), size(b,2))

      ! my vars      
      logical :: dSave = .True.
      integer :: info

      integer :: ipivot(size(b))
      real(8), allocatable :: A0(:,:), b0(:,:)
      
      ! check dims
      if (size(A,1) /= size(b,1)) then
         ! return. set info if required
         if(present(info_)) then
            info_ = -1
         end if
         return
      end if
      
      ! set dSave
      if (present(dSave_)) then
         dSave = dSave_
      end if
            
      ! save data if required
      if (dSave) then
         allocate( A0(size(A,1), size(A,2)) )
         allocate( b0(size(b,1), size(b,2)) )         
         A0 = A
         b0 = b
      end if
            
      ! int n = number of linear equations, i.e., order of matrix A
      ! int nrhs = number of right hand sides in b
      ! real a(n,n) =  in: coefficient matrix A.
      !             out: factors L and U from the factorization
      ! int lda = n leading dimension of a
      ! int ipv(n) = pivote vector for fact
      ! real b(n) = rhs
      ! int ldb = n  leading dimension of b
      ! int info =  exit code
      
      call dgesv(size(a,1), size(b,2), a, size(a,1), ipivot, b, size(b,1), info)                         
      
      ! return value
      x = b
      
      ! restore data if required
      if (dSave) then         
         A = A0
         b = b0
      end if
      
      ! return info if required
      if(present(info_)) then
         info_ = info
      end if
         
   end function linalg_solve_bn
#endif
      
	! ---
	   
#ifdef LAPACK
   !> Computes the solution to a real system of \(n\) linear equations \(A x = b\) (LAPACK)
   function linalg_solve_b1(A, b, dSave_, info_) result (x)
      real(8), intent (in out) :: a(:,:)  !! A(n,n) = coefficients matrix 
      real(8), intent (in out) :: b(:)    !! b(n) = RHS matrix, one column 
      logical, optional, intent(in)  :: dSave_  !! preserve original input data
      integer, optional, intent(out) :: info_   !! operation info
      real(8) :: x(size(b))
      
      ! my vars
      logical :: dSave = .True.
      integer :: info
      real(8) :: bk(size(b),1), xk(size(b),1)
      
      ! set dSave
      if (present(dSave_)) then
         dSave = dSave_
      end if
      
      ! set pars
      bk(:,1) = b
      
      ! process
      xk = linalg_solve_bn(A, bk, dSave, info)
      
      if (present(info_)) then
         info_ = info
      end if      
      
      ! return
      x = xk(:,1)
      
   end function linalg_solve_b1
#endif

	! ---
   
#ifdef FGSL
   !> 1-D interpolation using FGSL ([GNU Scientific Library]())
   function interp_interp1d(xi, x, y, itype_) result(yi)
   
      use fgsl
      
      ! arguments
      real(8) :: xi(:)        !! x locations for computing interpolations
      real(8) :: x(:), y(:)   !! known data
      real(8) :: yi(size(xi)) !! interpolated values
      character(*), optional :: itype_ 
                              !! interpolation type. default=`cspline`.
                              !! valid: `linear`, `poly`, `cspline`, `akima`, `steffen`
                              !! `cspline_p`, `akima_p`
      
      ! local vars
      character(20) :: itype
      real(8) :: yii(size(yi))
      type(fgsl_spline) :: aSpline
      type(fgsl_interp_accel) :: acc
      integer(fgsl_int) :: status
      integer(fgsl_size_t) :: n
      type(fgsl_interp_type) :: itypek
      integer :: i
      
      ! set itype
      if (present(itype_)) then
         itype = itype_
      end if
      
      select case (trim(itype))
      
         case('linear')
            itypek = fgsl_interp_linear
         case('poly')
            itypek = fgsl_interp_polynomial
         case('cspline')
            itypek = fgsl_interp_cspline
         case('akima')
            itypek = fgsl_interp_akima
         case('steffen')
            itypek = fgsl_interp_steffen
         case('cspline_p')
            itypek = fgsl_interp_cspline_periodic
         case('akima_p')
            itypek = fgsl_interp_akima_periodic
         case default
            print *, '--> you suck -> cspline'
            itypek = fgsl_interp_cspline
            
      end select
            
      n = size(x)    
      acc = fgsl_interp_accel_alloc()
      aSpline = fgsl_spline_alloc( itypek, n )
      status = fgsl_spline_init ( aSpline, x, y )
      
      ! interpolate xi     
      do i=1, size(yii)
         yii(i) = fgsl_spline_eval (aSpline, xi(i), acc)
      end do
      
      call fgsl_spline_free(aSpline)
      call fgsl_interp_accel_free(acc)
      
      ! return
      yi = yii

   end function interp_interp1d
#endif

   ! ---

#ifdef FFTW
   !> Computes fft 1d, real*8 -> complex*8, zero-padded N/1 -> N (FFTW-3).
   function fft_1d_r2c(y,n_) result(yfft)
      
      use iso_c_binding
      include 'fftw3.f03'
      
      ! args
      real(8) :: y(:)
      integer, optional, intent(in) :: n_
      complex(8), allocatable :: yfft(:)
      
      ! my vars
      integer(8) :: fPlan
      integer :: n
      
      ! set n
      n = size(y)
      if ( present(n_) ) n = n_
      
      ! allocate and set return
      allocate( yfft(n) )

      ! plan and execute fftw
      call dfftw_plan_dft_r2c_1d(fPlan, n, y, yfft, FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(fPlan, y, yfft)
      call dfftw_destroy_plan(fPlan)
   
   end function fft_1d_r2c
#endif

   ! ---

#ifdef FFTW
   ! fft (FFTW-3)
   ! computes ifft id, complex*8 -> real*8
   ! FFTW_PRESERVE_INPUT :
   !     - default, except for c2r and hc2r (causes worse performance)
   !     - for multi-dimensional c2r transforms, not available -> run error
   function ifft_1d_c2r(yFFT, dSave_) result(y)
      use iso_c_binding
      include 'fftw3.f03'
      
      ! args
      complex(8)  :: yFFT(:)
      logical, optional :: dSave_
      real(8), allocatable :: y(:)
      
      ! my vars
      complex(8), allocatable :: yFFT0(:)
      ! real(8), allocatable :: y0(:)
      integer(8) :: fPlan
      integer :: n
      logical :: dSave = .True.
      
      ! set dSave
      if (present(dSave_)) then
         dSave = dSave_
      end if
            
      ! set n
      n = size(yFFT)
      
      ! save data if required
      if (dSave) then
         allocate( yFFT0(n) )
         yFFT0 = yFFT
      end if
      
      allocate( y(n) ) 
      
      ! ifftw
      call dfftw_plan_dft_c2r_1d(fPlan, n, yFFT, y, FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(fPlan, yFFT, y)
      call dfftw_destroy_plan(fPlan)
      
      ! scale down
      y = y / n

      ! restore data if required
      if (dSave) then
         yFFT = yFFT0
      end if

      ! cleanup
      deallocate(yFFT0)
      
   end function ifft_1d_c2r
#endif

   ! ---

#ifdef NETCDF   
   !> Check I/O status of a netcdf operation
   subroutine nc_check(err)
      use netcdf
      integer, intent (in) :: err

      if(err /= nf90_noerr) then
         print *, '** NetCDF Error: ', nf90_strerror(err)
         stop "Stopped"
      end if
   end subroutine nc_check
#endif
   
#ifdef NETCDF         
   !> Get dimension lengths for the given NetCDF file and var 
   subroutine io_nc_getvardims(fname, vname, dlens, ndims_, dnames_, vid_, sort_)
   
      use netcdf
      implicit none
      
      ! parameters
      character(*), intent(in)            :: fname       !! NetCDF file name
      character(*), intent(in)            :: vname       !! Var name
      integer, allocatable, intent(out)   :: dlens(:)    !! Dimensions length
      integer, optional, intent(out)      :: ndims_      !! No. of dimensions
      character(*), allocatable, optional, intent(out) &
                                          :: dnames_(:)  !! Dimensions name
      integer, optional, intent(out)      :: vid_        !! Var id
      logical, optional, intent(in)       :: sort_       !! Sort dimensions to e.g. (..,t,lat,lon)
      
      ! local vars
      integer :: ncid                  ! nc file id
      integer :: vid, vndims0 ! nc var id, var ndims
      character(nf90_max_name) :: atxt ! any text  
      integer, dimension(nf90_max_var_dims) :: vdimids
      integer :: i
      logical :: sortdims = .true.
   
      ! open file
      call nc_check( nf90_open(trim(fname), nf90_nowrite, ncid) )
   
      ! get var id
      call nc_check( nf90_inq_varid(ncid, trim(vname), vid) )
      if (present(vid_)) vid_ = vid

      ! get no of dimensions and ids
      call nc_check( nf90_inquire_variable(ncid, vid, ndims=vndims0, dimids=vdimids) )
      
      ! set ndims if present
      if (present(ndims_)) then
         ndims_ = vndims0
      end if
      
      ! set dnames if required
      if (present(dnames_)) allocate(dnames_(vndims0))
      
      ! allocate dlens
      allocate(dlens(vndims0))

      ! get dims lens x-y, (it is not row-cols!)
      do i=1, vndims0
         call nc_check( nf90_inquire_dimension(ncid, vdimids(i), name=atxt, len=dlens(i)) )
         if (present(dnames_)) dnames_(i) = trim(atxt)
      end do

      ! arrange dims info from internal (dim1..dimn) to regular order (dimn..dim1) order
      ! e.g. (lon,lat,t) -> (t,lat,lon)
      ! e.g. (180, 170, 24)  -> (24, 170, 180)
      ! sys_ = internal system usage (called from internal level routines)
      if (present(sort_)) sortdims = sort_

      if (sortdims) then
         call numfor_flip(dlens)
         if (present(dnames_)) call numfor_flip(dnames_)
      end if
   
      ! bye
      call nc_check( nf90_close(ncid) )
   

   end subroutine io_nc_getvardims
#endif
      
   ! ---

#ifdef NETCDF         
   !> Reads a 1-D, `real(8)` var from a NetCDF file.
   ! file is open and closed transparently.
   subroutine io_nc_getvar_1d_r8(fname, vname, var, start_, count_, err_)

      use netcdf
      implicit none

      ! parameters
      character(*), intent(in) :: fname   !! NetCDF file name
      character(*), intent(in) :: vname   !! variable name to read
      real(8), allocatable, intent(out) :: var(:)  !! array save var data
      integer, optional, intent(in) :: start_      !! data reading start
      integer, optional, intent(in) :: count_      !! data count
      integer, optional, intent(out) :: err_       !! error info

      ! local vars
      integer :: ndims, ncid, varid
      integer :: start0, count0
      integer, allocatable :: dimsl(:)

      ! main

      ! get var dims
      call io_nc_getvardims(fname, vname, dimsl, ndims_=ndims, vid_=varid)

      ! set defaults
      start0 = 1
      count0 = dimsl(1)
 
      ! set start from params
      if (present(start_)) then
         start0  = start_
         count0 =  dimsl(1) - start0 + 1 
      end if

      ! set count from params
      if (present(count_)) count0 = count_

      ! allocate mem
      allocate(var(count0))

      ! open file
      call nc_check( nf90_open(fname, nf90_nowrite, ncid) )

      ! locate var
      call nc_check( nf90_inq_varid(ncid, vname, varid) )

      ! get data
      call nc_check( nf90_get_var(ncid, varid, var, start=[start0], count=[count0]) )

   end subroutine io_nc_getvar_1d_r8
#endif

   ! ---

#ifdef NETCDF         
   ! reads a 2-D, real(4) var from fname, using star_ and count_ if provided
   ! requires: fname, vname, var
   ! returns: var, err_ if provided
   ! file is opened and closed transparently
   subroutine io_nc_getvar_2d_r4(fname, vname, var, start_, count_, err_)

      use netcdf
      implicit none

      ! parameters
      character(*), intent(in) :: fname
      character(*), intent(in) :: vname
      real(4), allocatable, intent(out) :: var(:,:)
      integer, optional, intent(in)     :: start_(:)
      integer, optional, intent(in)     :: count_(:)
      integer, optional, intent(out)    :: err_

      ! local vars
      real(4), allocatable :: var0(:,:)
      integer :: ndims, ncid, varid, j
      integer, allocatable :: start0(:), count0(:), vdimsl(:)

      ! main
      print *, '  -- io_nc_getvar_2d_r4 : ', trim(vname)

      ! get var dims
      call io_nc_getvardims(fname, vname, vdimsl, ndims_=ndims, vid_=varid, sort_=.false.)
      print *, '     vdimsl = ', vdimsl

      ! set default start 
      allocate(start0(ndims))
      start0 = 1

       ! set start from params
      if (present(start_)) then
         start0(1:size(start_)) = start_
         call numfor_flip(start0)
      end if

      ! set default count
      ! by default, read all data
      ! another option: count(:numDims) = shape(values) and count(numDims + 1:) = 1
      allocate(count0(ndims))
      count0 = vdimsl

      ! set count from params
      if (present(count_)) then
         count0 = count_
         call numfor_flip(count0)
      end if

      ! allocate internal 2d var, (dim1 ... dimn)
      allocate( var0(count0(1),count0(2)) )

      ! open file
      call nc_check( nf90_open(fname, nf90_nowrite, ncid) )

      ! locate var
      call nc_check( nf90_inq_varid(ncid, vname, varid) )

      ! get data
      print *, '     getting data, shape =', shape(var0)
      call nc_check( nf90_get_var(ncid, varid, var0, start=start0, count=count0) )

      ! allocate final 2d var, (dimn .. dim1)
      ! final var contains data flipped from dim1..dimn to dim2..dim1
      allocate( var(count0(2),count0(1)) )

      print *, '     assembling to shape ', shape(var0)
      
      ! arrange output ordered from dim1.. dim3 to dim3..dim1
      do j=1, count0(2)
         var(j,:) = var0(:,j)
      end do

      ! deallocate var0
      deallocate(var0)

      ! close file
      call nc_check( nf90_close(ncid) )

   end subroutine io_nc_getvar_2d_r4
#endif

   ! ---

#ifdef NETCDF            
   !> Reads a 3-D, real(4) var from `fname`, using `start_` and `count_` if provided.
   ! File is opened and closed transparently.
   subroutine io_nc_getvar_3d_r4(fname, vname, var, start_, count_, err_)

      use netcdf
      implicit none

      ! parameters
      character(*), intent(in) :: fname   !! NetCDF file name
      character(*), intent(in) :: vname   !! Variable name to read
      real(4), allocatable, intent(out) :: var(:,:,:) !! Variable
      integer, optional, intent(in)     :: start_(:)  !! Start of file reading
      integer, optional, intent(in)     :: count_(:)  !! No. of elements to read
      integer, optional, intent(out)    :: err_       !! Error code (not implemented yet)

      ! local vars
      real(4), allocatable :: var0(:,:,:)
      integer :: ndims, ncid, varid, j, k
      integer, allocatable :: start(:), count(:), vdimsl(:)

      ! main

      ! get var dims
      call io_nc_getvardims(fname, vname, vdimsl, ndims_=ndims, vid_=varid, sort_=.false.)

      ! set default start 
      allocate(start(ndims))
      start = 1

       ! set start from params
      if (present(start_)) then
         start(1:size(start_)) = start_
         call numfor_flip(start)
      end if

      ! set default count
      ! by default, read all data
      ! another option: count(:numDims) = shape(values) and count(numDims + 1:) = 1
      allocate(count(ndims))
      count = vdimsl

      ! set count from params
      if (present(count_)) then
         count = count_
         call numfor_flip(count)
      end if

      ! allocate internal 3d var, (dim1 ... dimn)
      allocate( var0(count(1),count(2),count(3)) )

      ! open file
      call nc_check( nf90_open(fname, nf90_nowrite, ncid) )

      ! locate var
      call nc_check( nf90_inq_varid(ncid, vname, varid) )

      ! get data

      call nc_check( nf90_get_var(ncid, varid, var0, start=start, count=count) )

      ! allocate final 3d var, (dimn .. dim1)
      ! final var contains data flipped from dim1..dimn to dim2..dim1
      allocate( var(count(3),count(2),count(1)) )
      
      ! arrange output ordered from dim1.. dim3 to dim3..dim1
      do k=1, count(3)
         do j=1, count(2)
            var(k,j,:) = var0(:,j,k)
         end do
      end do

      ! deallocate var0
      deallocate(var0)

      ! close file
      call nc_check( nf90_close(ncid) )

   end subroutine io_nc_getvar_3d_r4
#endif

	! ---
	
#ifdef NETCDF            	
   !> Gets the value of attribute `attname` (1-D, real(4))
   subroutine io_nc_getatt_1d_r4(fname, vname, attname, attval)

      use netcdf
      implicit none

      ! parameters
      character(*), intent(in) :: fname      !! NetCDF file name
      character(*), intent(in) :: vname      !! Var name
      character(*), intent(in) :: attname    !! Att name
      real(4), allocatable, intent(out) :: attval(:)  !! Att value

      ! local vars
      integer :: ncid, varid, attype, attlen, status

      ! open file
      call nc_check( nf90_open(fname, nf90_nowrite, ncid) )

      ! get var id
      call nc_check( nf90_inq_varid(ncid, vname, varid) )

      ! get att value size
      call nc_check( nf90_inquire_attribute(ncid, varid, attname, attype, attlen) )

      ! check data type
      if (attype /= NF90_FLOAT) then
         call nc_check(NF90_EBADTYPE)
      end if

      ! allocate space to hold attribute values
      allocate(attval(attlen), stat=status)

      if (status /= 0 ) then
         stop "** Error: memory could not be allocated for attributes"
      end if

      ! Read the attributes.
      call nc_check( nf90_get_att(ncid, varid, attname, attval) )

      ! close file
      call nc_check( nf90_close(ncid) )

   end subroutine io_nc_getatt_1d_r4
#endif

	! ---

#ifdef NETCDF            	
     !> Gets the value of attribute `attname` (escalar, real(4))
   subroutine io_nc_getatt_r4(fname, vname, attname, attval)

      use netcdf
      implicit none

      ! parameters
      character(*), intent(in) :: fname      !! NetCDF file name
      character(*), intent(in) :: vname      !! Var name
      character(*), intent(in) :: attname    !! Att name
      real(4), intent(out)     :: attval     !! Att value

      ! local vars
      real(4), allocatable :: attval0(:)

      call io_nc_getatt_1d_r4(fname, vname, attname, attval0)

      ! set out value
      attval = attval0(1)

   end subroutine io_nc_getatt_r4
#endif


end module
