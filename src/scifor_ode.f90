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
        real(8), allocatable, public :: t(:), y(:,:), troot(:), yroot(:,:)
        integer, allocatable, public  :: nfroot(:,:)
        real(8), allocatable  :: rtol(:), atol(:), rwork(:)
        integer               :: neq, itol, itask, istate, iopt, jtype, jac, nfr
        integer, allocatable  :: iwork(:)                
        character(20)         :: solver
        logical               :: initialized = .False.
        integer               :: nchunk 
        
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
   ! To provide a canonical constructor function to be called as `ode = Odesolv()`
   ! However, it is not called intrinsecally when the object is constructed
   ! It is up to the user to call it
   ! To work around, the function `Odesolv%init()` is provided and is called by the
   ! constructor as only task. This way, `init()` can be called when necessary, as 
   ! when the `solve()` procedure is called.
   ! Constructor kept only for style completness
   !
   function constructor() result(Ode)
      type(Odesolv)	:: Ode
            
      ! init
      call Ode%init()

   end function constructor


   !> Init
   !
   subroutine init(this, neq)
      class(Odesolv)    :: this
      integer, optional :: neq
      integer           :: neqw

      ! defaults
      neqw = 1
      if (present(neq)) neqw = neq

      ! init values
      this%solver = 'sodar'   ! default solver
      this%neq    = neqw      ! n of eq to solve
      this%itol   = 1         ! rtol, atol scalars (not arrays)
      this%itask  = 5         ! auto time step dt
      this%istate = 1         ! first call
      this%iopt   = 1         ! optional inputs (rwork, iwork) present (tcritic)
      this%nfr = 0            !  n of funct roots
      this%jtype  = 2         ! computed jacobian
      this%nchunk = 1024      ! n elem to alloc per chunk (t, y, roots, )

      ! abs error tolerance. default: scalar, no rtol
      if (allocated(this%atol)) deallocate(this%atol)
      this%atol = [1d-6] 

      ! rel error tolerance (scalar)
      if (allocated(this%rtol)) deallocate(this%rtol)
      this%rtol = [0d0]

      ! rwork
      if (allocated(this%rwork)) deallocate(this%rwork)
      allocate(this%rwork(this%neq*128))    ! experimental value !!!
      this%rwork  = 0d0

      ! iwork
      if (allocated(this%iwork)) deallocate(this%iwork)
      allocate(this%iwork(this%neq+25))
      this%iwork  = 0

      ! clean results
      if (allocated(this%t)) deallocate(this%t)
      if (allocated(this%y)) deallocate(this%y)
      if (allocated(this%troot)) deallocate(this%troot)
      if (allocated(this%yroot)) deallocate(this%yroot)
      if (allocated(this%nfroot)) deallocate(this%nfroot)

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
      procedure()           :: fdy     ! diff function
      real(8)               :: y0(:)   ! initial cond
      real(8), optional     :: trng(:), tpts(:) ! t range
      real(8), optional     :: rtol(:), atol(:) ! t points
      procedure(), optional :: jac, froots      ! jac() and froots() functions
      integer, optional     :: nfroots              ! n function rules in froots()
      real(8), optional     :: dt               ! desired dt step

      ! local vars
      real(8), allocatable, target  :: t1(:), t2(:), y1(:,:), y2(:,:)
      real(8),  pointer             :: pt(:), py(:,:)
      real(8), allocatable, target  :: troot1(:), troot2(:), yroot1(:,:), yroot2(:,:)
      real(8), pointer              :: ptroot(:), pyroot(:,:)
      integer, allocatable, target  :: nfroot1(:,:), nfroot2(:,:)
      integer, pointer              :: pnfroot(:,:)
      !
      real(8)     :: tw, tout, yw(size(y0)), dtw
      integer, allocatable :: nfrootw(:)
      integer     :: it, ir
      logical     :: doAnother
            
      ! ** number of functions is defined using the number of initial conditions `y0`**

      ! validate t
      if (.not. present(trng) .and. .not. present(tpts)) then
         stop '** Error: Please provide either time range `trng` or time points `tpts`'
      end if
      !
      if (present(trng) .and. present(tpts)) then
         stop '** Error: Please provide either time range or time points only'
      end if

      ! validate froots
      if (present(froots) .and. .not. present(nfroots)) then
         stop '** Error: Please provide n of roots `nfr` if `froots()` is provided'
      end if

      ! validate dt (for trange, dt=auto)
      dtw = 1d-3     
      if (present(dt)) dtw = dt


      ! ok, init object
      call this%init(size(y0))
      
      ! 1st call y = y0 (initial conditions)
      yw = y0

      ! if time range provided, set up
      if (present(trng)) then
      
         ! check size
         if (size(trng) /= 2) then
            print *, '** Error: Plese provide time range `trng` as a 2-value array'
            stop
         end if

         ! 1st call t = trng(1)
         tw = trng(1)         

         ! set tstop
         this%rwork(1) = trng(2)            
         
         ! 1st desired tout
         tout = tw + dtw
      end if

      ! if time points provided, set up
      if (present(tpts)) then
      
         ! check size
         if (size(tpts) < 2) then
            print *, '** Error: Plese provide time points `tpts` at least as a 2-value array'
            stop
         end if
         
         this%itask = 1  ! manual dt: solve odt at t=tpts(i)
         this%iopt  = 0  ! no options present in rwork, iwork  
         
         ! 1st call
         tw         = tpts(1)
         tout       = tpts(2)
      end if

      ! use jac if provided (normal jacobian)
      if (present(jac)) then
         this%jtype = 1
      end if

      ! set rtol
      if (present(rtol)) then
         if (size(rtol) /= 1 .and. size(rtol) /= this%neq) then
            stop "** Error: relative tolerance `rtol` lenght should be 1 or `neq`"
         end if
         
         if (size(this%rtol) /= size(rtol)) deallocate(this%rtol)
         this%rtol = rtol
      end if

      ! set atol
      if (present(atol)) then
         if (size(atol) /= 1 .and. size(atol) /= this%neq) then
            stop "** Error: absolute tolerance `atol` lenght should be 1 or `neq`"
         end if
         
         if (size(this%atol) /= size(atol)) deallocate(this%atol)
         this%atol = atol
      end if

      ! set itol flag
      if ( size(this%rtol) == 1 .and. size(this%atol) == 1) then
         this%itol = 1
      else if (size(this%rtol) == 1 .and. size(this%atol) > 1) then
         this%itol = 2
      else if (size(this%rtol) > 1 .and. size(this%atol) == 1) then
         this%itol = 3
      else
         this%itol = 4
      end if
     
      ! allocate nfroot if froots() provided
      if (present(froots)) then
         this%nfr = nfroots
         allocate(nfrootw(nfroots))
      end if

      ! initial allocation for results arrays (dynamic growing @ nchunk elements!)
      allocate (t1(this%nchunk))
      allocate (y1(this%nchunk, this%neq))
      
      ! init pointers
      pt => t1
      py => y1

      ! same if froots() provided
      if (present(froots)) then
         allocate (troot1(this%nchunk))
         allocate (yroot1(this%nchunk, this%neq))
         allocate (nfroot1(this%nchunk,this%nfr))

         ptroot => troot1
         pyroot => yroot1
         pnfroot => nfroot1
      end if
  
      ! init loop, one step per dt 
      
      it  = 0
      ir = 0
      doAnother = .True.

      do while(doAnother)
      
         call dlsodar( fdy, size(y0), yw, tw, tout,         &
                     this%itol, this%rtol, this%atol,       &
                     this%itask, this%istate, this%iopt,    &
                     this%rwork, size(this%rwork),          &
                     this%iwork, size(this%iwork),          &
                     jac, this%jtype,                       &
                     froots, this%nfr, nfrootw )

         ! print *, it, tout, tw, yw

         if (this%istate < 0) then
            print *, '** Error: ISTATE =', this%istate
            stop
         end if


         select case (this%istate)
         
            ! normal result, no root: advance t
            case (2)
               it = it+1
            
               ! check pt/py arrays occupation; increase if required
               if (it > size(pt)) then
                  if (associated(pt, target=t1)) then

                     ! allocate space in t2/y2
                     allocate (t2(it-1 + this%nchunk))
                     allocate (y2(it-1 + this%nchunk, this%neq))

                     ! copy t1,y1 -> t2,y2
                     t2(:it-1) = t1 
                     y2(:it-1,:) = y1
                     
                     ! update pointers
                     pt => t2
                     py => y2
                     
                     ! release t1,y1
                     deallocate(t1)
                     deallocate(y1)

                  else ! pt2 is allocated

                     allocate (t1(it-1 + this%nchunk))
                     allocate (y1(it-1 + this%nchunk, this%neq))

                     t1(:it-1) = t2 
                     y1(:it-1,:) = y2

                     pt => t1
                     py => y1
                     
                     deallocate(t2)
                     deallocate(y2)

                  end if
               end if

               ! save data		
               pt(it) = tw
               py(it,:) = yw

               ! another loop?
               select case (this%itask)
            
                  ! t = fixed points, tout=set, dt=manual
                  case (1) 
                     if (it < size(tpts)-1) then
                        tout = tpts(2+it)
                     else
                        doAnother = .False.
                     end if
            
                  ! t = range, tout=auto, dt=auto
                  case (5) 
                     doAnother = tw < trng(2)
               
               end select

                     
            ! special result, root found
            case (3)
               ir = ir + 1
               ! print *, '# -- root found :', tw, ir, nfrootw

               ! check ptroot/pyroot/pnfroot arrays occupation; increase if required
               if (ir > size(ptroot)) then
                  if (associated(ptroot, target=troot1)) then

                     ! allocate space in troot2/yroot2/nfroot2
                     allocate (troot2(ir-1 + this%nchunk))
                     allocate (yroot2(ir-1 + this%nchunk, this%neq))
                     allocate (nfroot2(ir-1 + this%nchunk, this%nfr))
                     
                     ! copy troot1,yroot1,nfroot1 -> troot2,yroot2,nfroot2
                     troot2(:ir-1)    = troot1 
                     yroot2(:ir-1,:)  = yroot1
                     nfroot2(:ir-1,:) = nfroot1

                     ! update pointers
                     ptroot  => troot2
                     pyroot  => yroot2
                     pnfroot => nfroot2
                     
                     ! release mem
                     deallocate(troot1)
                     deallocate(yroot1)
                     deallocate(nfroot1)

                  else ! troot2 is allocated

                     ! allocate space in troot1/yroot1/nfroot1
                     allocate (troot1(ir-1 + this%nchunk))
                     allocate (yroot1(ir-1 + this%nchunk, this%neq))
                     allocate (nfroot1(ir-1 + this%nchunk, this%nfr))
                     
                     ! copy troot2,yroot2,nfroot2 -> troot1,yroot1,nfroot1
                     troot1(:ir-1)    = troot2
                     yroot1(:ir-1,:)  = yroot2
                     nfroot1(:ir-1,:) = nfroot2

                     ! update pointers
                     ptroot  => troot1
                     pyroot  => yroot1
                     pnfroot => nfroot1
                     
                     ! release mem
                     deallocate(troot2)
                     deallocate(yroot2)
                     deallocate(nfroot2)

                  end if
               end if

               ! save data		
               ptroot(ir) = tw			
               pyroot(ir,:) = yw
               pnfroot(ir,:) = nfrootw
            
         end select
         
         ! set up next iter
         this%istate = 2
  
      end do
      
      ! save results
      ! a: t, u
      this%t = pt(1:it)
      this%y = py(1:it,:)
      !
      ! b: roots
      if (present(froots)) then
         this%troot = ptroot(1:ir)
         this%yroot = pyroot(1:ir,:)
         this%nfroot = pnfroot(1:ir,:)
      end if

      ! release mem
      if (allocated(t1)) deallocate(t1)
      if (allocated(t2)) deallocate(t2)
      if (allocated(y1)) deallocate(y1)
      if (allocated(y2)) deallocate(y2)
      !     
      if (present(froots)) then
         if (allocated(troot1)) deallocate(troot1)
         if (allocated(troot2)) deallocate(troot2)
         if (allocated(yroot1)) deallocate(yroot1)
         if (allocated(yroot2)) deallocate(yroot2)
         if (allocated(nfroot1)) deallocate(nfroot1)
         if (allocated(nfroot2)) deallocate(nfroot2)
      end if

      ! print *, ''
      ! print *, '-- main loop done'
      ! print *, ''

      ! print *, ' --- t: shape=', shape(this%t), ', it=', it
      ! print *, this%t
      ! print *, ''
      
      ! print *, ' --- y: shape=', shape(this%y), ', it=', it
      ! call print_mat_r8(this%y)
      ! print *, ''
            
      ! if (present(froots)) then
      !    print *, ''
      !    print *, ' +-- troot: size=', size(this%troot), ', ir=', ir
      !    print *, this%troot
      !    print *, ''

      !    print *, ' --- yroot: size=', shape(this%yroot), ', ir=', ir
      !    call print_mat_r8( this%yroot )
      !    print *, ''

      !    print *, ' --- nroot: size=', shape(this%nfroot), ', ir=', ir
      !    call print_mat_r8( 1d0*this%nfroot, 'F12.0' )
      !    print *, ''
      ! end if

    end subroutine solve
    
    
    subroutine print_mat_r8(A, fmt_)
      
      real(8), intent(in) :: A(:,:)                !! matrix to print
      character(*),optional, intent(in) :: fmt_    !! print format (default:'ES12.4').
      
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


end module scifor_ode