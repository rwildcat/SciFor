!> #A modern Fortran interface to ODE libraries
!
! ##Initial Value Problem (IVP)
!
! ###[ODEPACK](https://www.netlib.org/odepack/opkd-sum),
! a Systematized Collection of ODE Solvers.
!
! ODEPACK is a collection of Fortran solvers for the initial value
! problem for ordinary differential equation systems.  It consists of nine
! solvers, namely a basic solver called `LSODE` and eight variants of it --
! LSODES, LSODA, LSODAR, LSODPK, LSODKR, LSODI, LSOIBT, and LSODIS.
! The collection is suitable for both stiff and nonstiff systems.  It
! includes solvers for systems given in explicit form, dy/dt = f(t,y),
! and also solvers for systems given in linearly implicit form, 
! A(t,y) dy/dt = g(t,y).  Two of the solvers use general sparse matrix
! solvers for the linear systems that arise.  Two others use iterative
! (preconditioned Krylov) methods instead of direct methods for these
! linear systems.  The most recent addition is LSODIS, which solves
!
! ## Boundary Value Problem (BVP) Algorithms -- ** To Be Implemented **
!
! * See `L.F. Shampine, P. H. Muir and H. Xu, “A User-Friendly Fortran BVP Solver”`.
!
!     This algorithm may be used by [`scipy.integrate.solve_bvp`](https://bit.ly/2WQeQOH)
!
! To compile during tests
!       $ gfortran ../src/scifor_ode.f90 ode_scifor-test.f90
!
module scifor_ode

    implicit none
    private

    type, public :: Odesolv
        private
        integer :: neq, itol, itask, istate, iopt, jtype, jac, nroots
        real(8), allocatable :: rtol(:), atol(:), rwork(:)
        real(8), allocatable :: dt(:), dy(:,:)
        integer, allocatable :: iwork(:), jroot(:,:)
        character(20)        :: solver
        logical              :: initialized = .False.
        integer              :: nchunk 
        
    contains
        private
        procedure, public :: init
        procedure, public :: print
        procedure, public :: solve

    end type Odesolv

    ! constructor
    interface Odesolv
        module procedure constructor
    end interface

contains

   !> Constructor
   function constructor(neq) result(Ode)

      ! params
      type(Odesolv)     :: Ode
      integer,optional  :: neq

      ! local vars
      integer           :: neq_ = 1

      if (present(neq)) neq_ = neq
      call Ode%init(1)

   end function constructor


   !> Init for normal, simpler oper
   subroutine init(this, neq)

      ! params
      class(Odesolv)    :: this
      integer, optional :: neq

      ! local vars
      integer           :: neq_ = 1

      ! chack neq
      if (present(neq)) neq_ = neq

      ! init
      this%solver = 'sodar'   ! default solver
      this%neq    = neq_         ! 1 eq
      this%itol   = 1         ! rtol=scalar, atol=scalar
      this%itask  = 5         ! auto time step
      this%istate = 1         ! first call
      this%iopt   = 1         ! optional inputs (rwork, iwork) present (tcritic)
      this%nroots = 0         ! no roots required
      this%jtype  = 2         ! internal jacobian
      this%nchunk = 1024      ! n elem to alloc

      ! neq
      if (present (neq)) neq_ = neq

      this%neq = neq_

      ! abs error tolerance (scalar)
      if (allocated(this%atol)) deallocate(this%atol)
      allocate(this%atol(1))
      this%atol = 1d-6  

      ! rel error tolerance (scalar)
      if (allocated(this%rtol)) deallocate(this%rtol)
      allocate(this%rtol(1))
      this%rtol = 0d0

      if (allocated(this%rwork)) deallocate(this%rwork)
      allocate(this%rwork(this%neq*128))    ! experimental !!!
      this%rwork  = 0d0

      if (allocated(this%iwork)) deallocate(this%iwork)
      allocate(this%iwork(this%neq+25))
      this%iwork  = 0

      ! system roots
      if (allocated(this%jroot)) deallocate(this%jroot)

      ! release dt, dy
      if (allocated(this%dt)) deallocate(this%dt)
      if (allocated(this%dy)) deallocate(this%dy)

      this%initialized = .True.
      
   end subroutine init


   !> Just print some info
   subroutine print(this)
      class(Odesolv) :: this
      print *, this%initialized, this%neq, this%itol, this%rtol, this%atol, this%jtype, trim(this%solver)
   end subroutine print


   !> Solve ODE (LSODAR)
   subroutine solve(this, fdy, y0, trng, tpts, rtol, atol, jac, froots, nfroots, dt)
      class(Odesolv)        :: this
      procedure()           :: fdy
      real(8)               :: y0(:)
      real(8), optional     :: trng(:), tpts(:), rtol(:), atol(:)
      procedure(), optional :: jac, froots
      integer, optional     :: nfroots
      real(8), optional     :: dt

      ! local (working) vars
      procedure(), pointer :: froots0
      real(8), allocatable :: dtw(:), dyw(:,:)
      real(8)     :: yw(size(y0)), tw, dt0, tout
      integer     :: itot, ichnk, iroot
      integer, allocatable :: jrootw(:)
      logical     :: doAnother

      ! validate t
      if (.not. present(trng) .and. .not. present(tpts)) then
         stop '** Error: Please provide either time range `trng` or time points `tpts`'
      end if
      !
      if (present(trng) .and. present(tpts)) then
         stop '** Error: Please provide only either time range or time points only'
      end if

      ! validate froots
      if (present(froots) .and. .not. present(nfroots)) then
         stop '** Error: Please provide `nroots` if `froots` provided'
      end if

      ! init object (default: auto step, mode 5)
      call this%init(size(y0))
      yw = y0

      ! if auto dt mode set tcritic 
      if (present(trng)) then
         ! check size
         if (size(trng) /= 2) then
            print *, '** Error: Plese provide time range `trng` as a 2-value array'
            stop
         end if

         tw = trng(1)         
         this%rwork(1) = trng(2)            

         dt0 = 1d-3      ! default initial dt for auto mode
         if (present(dt)) dt0 = dt
         
         tout = tw + dt0     
      end if

      ! if manual dt, set manual mode
      if (present(tpts)) then
         this%itask = 1  ! manual dt: solve odt(tout)
         this%iopt  = 0  ! no options
         tw         = tpts(1)
         tout       = tpts(2)
      end if

      ! use jac if provided
      if (present(jac)) then
         ! user-provided, normal (full) jacobian
         this%jtype = 1
      end if

      ! validate rtol
      if (present(rtol)) then
         if (size(rtol) /= size(this%rtol)) then
            deallocate(this%rtol)
            allocate(this%rtol(size(rtol)))
         end if
         this%rtol = rtol
      end if

      ! validate atol
      if (present(atol)) then
         if (size(atol) /= size(this%atol)) then
            deallocate(this%atol)
            allocate(this%atol(size(atol)))
         end if
         this%atol = atol
      end if

      ! validate tol sizes
      if ( size(this%rtol) /= 1 .and. size(this%rtol) /= this%neq .or.  &
           size(this%atol) /= 1 .and. size(this%atol) /= this%neq ) then
            stop "** Error: Please provide either one (scalar) or full (array) `rtol`|`atol` values"
      end if

      ! set itol (default = 1)
      if ( size(this%rtol) == 1 .and. size(this%atol) == 1) then
         this%itol = 1
      else if (size(this%rtol) == 1 .and. size(this%atol) > 1) then
         this%itol = 2
      else if (size(this%rtol) > 1 .and. size(this%atol) == 1) then
         this%itol = 3
      else
         this%itol = 4
      end if
     
      ! allocate working arrays (dynamic growing)
      allocate(dtw(this%nchunk))
      allocate(dyw(this%nchunk, this%neq))
      allocate(this%jroot(this%nchunk,this%neq))

      ! allocate final and working jroot arrays
      if (present(froots)) then
         this%nroots = nfroots
         allocate(jrootw(nfroots))
      else
         allocate(jrootw(1))
      end if

      ! init loop
      itot  = 1
      ichnk = 1
      iroot = 1
      doAnother = .True.

      do while(doAnother)

         call dlsodar( fdy, size(y0), yw, tw, tout,         &
                     this%itol, this%rtol, this%atol,       &
                     this%itask, this%istate, this%iopt,    &
                     this%rwork, size(this%rwork),          &
                     this%iwork, size(this%iwork),          &
                     jac, this%jtype,                       &
                     froots, this%nroots, jrootw )

         print *, itot, tw, yw
  
         if (this%istate < 0) then
            print *, '** Error: ISTATE =', this%istate
            stop
         end if

         ! root found?
         if (this%istate .eq. 3) print *, '# -- root found :',  jrootw

         ! set up next iter
         this%istate = 2
         ichnk = ichnk + 1
         itot = itot + 1

         ! another loop?
         select case (this%itask)
            case (1)
               if (itot < size(tpts)) then
                  tout = tpts(itot+1)
               else
                  doAnother = .False.
               end if

            case (5)
               doAnother = tw < trng(2)
         end select
  
      end do


    end subroutine solve


end module scifor_ode