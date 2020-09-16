!> ode_scifor-test.f90
! Purpose    : Demonstrations of the SciFor interace to ODE solvers
! Author     : Ramon Solano
!              ramon.solano at gmail.com
!              Colima, Mexico
! To compile :
!  `$ gfortran ode_scifor-test.f90 -logpf `pkg-config --cflags --libs scifor`
!
program ode_scifor_test

   !use scifor_ode
   use scifor
   use ogpf
   implicit none
   
   type(Odesolv)  :: ode 
   type(Gpf)      :: plt
   
	! ode = Odesolv()

   ! ---

   ! example 1: flame radius
   !  1st deg, stiff + roots
   !  y' = y^2 - y^3
   !  y(0) = 2d-3
   !  t    = 0, 2/2d-3

   ! plot init
   call plt%title("Flame radious (stiff) - y' = y^2 - y^3")
   call plt%xlabel("t (s)")
   call plt%ylabel("r (m)")

   call ode%solve ( fdy, [2d-3], trng=[0d0, 2d0/2d-3] )
   call plt%plot(ode%t, ode%y(:,1))

   ! same, but also find roots according fgex()
   call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3], froots=fgex, nfroots=1 )
   call plt%plot(ode%t, ode%y(:,1), 'title "r"','', ode%troot, ode%yroot(:,1), 'w circles title "root at 5e-2"')

   ! same, but also t = specific points
   call ode%solve( fdy, [2d-3], tpts=[0d0, 1d1, 3d2, 4.5d2, 4.75d2, 4.9d2, 5d2, 5.1d2, 6d2, 10d2] )
   call plt%plot(ode%t, ode%y(:,1))

   ! same, but also find roots according to fgex()
   call ode%solve( fdy, [2d-3], tpts=[0d0, 1d1, 3d2, 4.5d2, 4.75d2, 4.9d2, 5d2, 5.1d2, 6d2, 10d2], froots=fgex, nfroots=1 )
   call plt%plot(ode%t, ode%y(:,1), 'title "y"','', ode%troot, ode%yroot(:,1), 'w circles title "root at 5e-2"')

   ! ---

   ! example 2: van der Pol oscillator
   !  2nd deg, stiff/no stiff: 
   !  y'' = µ(1 - y^2)y' - y
   !  y(0), y'(0) = [2, 0]
   !  t = [0, 100]
   
   ! plot
   call plt%title("van der Pol oscillator (stiff/nostiff) - y'' = µ(1 - y^2)y' - y")
   call plt%xlabel("x (m)")
   call plt%ylabel("y (m)")

   ! normal
   call ode%solve( fvdp, [2d0,0d0], trng=[0d0, 100d0] )
   call plt%plot(ode%t, ode%y(:,1))

   ! with roots
   call ode%solve( fvdp, [2d0,0d0], trng=[0d0, 100d0], froots=frvdp, nfroots=1 )
   call plt%plot(ode%t, ode%y(:,1), 'title "y"','', ode%troot, ode%yroot(:,1), 'w circles title "roots"')
   
   ! ---

   ! example 3: 3 chem substances ratios
   !  1st degree 3-eq system + roots + jacobian
   !  fchem3 = see below

   call ode%solve(                &
      fchem3, [1d0, 0d0, 0d0],    &
      tpts=[0d0, 4d-1, 4d0, 4d1, 4d2, 4d3, 4d4, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10],  &
      rtol=[1d-4], atol=[1.d-6, 1.d-10, 1.d-6],   &
      jac=fjchem3, froots=frchem3, nfroots=2 )

   ! show roots info
   print *, '3 chem substances ratio'
   print *, 'troot =', ode%troot
   print *, 'yroot ='
   call io_print_mat(ode%yroot)

   ! plot
   call plt%title("3 substances concentration ratio")
   call plt%xlabel("t (s)")
   call plt%ylabel("R (%)")
   call plt%loglog( &
      ode%t, ode%y(:,1),'title "S1"', '', &
      ode%t, ode%y(:,2),'title "S2"', '', &
      ode%t, ode%y(:,3),'title "S3"')

   
contains

   !> 1st degree, stiff ODE system - Flame radius 
   !		y' = y^2 - y^3
   ! subject to: 
   !		y(0) = delta
   !		0 <= t <= 2/delta
   subroutine fdy(neq,t,u,du)
      implicit none
      integer, intent(in)  ::  neq
      real(8), intent(in)  ::  t, u(neq)
      real(8), intent(out) :: du(neq)
      ! local vars
      real(8) :: y, dy(neq)
      !
      ! unpack primitives y, y', y", ...
      y  = u(1)
      ! compute desired derivatives y', y", ...
      dy(1) = y**2 - y**3
      ! assemble derivatives du=(y', y", ...)
      du = dy  
   end subroutine fdy
   !
   ! functions roots
   subroutine fgex (neq, t, y, ng, gout)
      integer :: neq, ng
      real(8) :: t, y(1), gout(1)
      gout(1) = y(1) - 5d-2
   end subroutine


   !> 2nd degree, no/yes stiff ODE system - van der Pol oscillator
   ! No stiff: mu small, e.g. mu = 1, t = [0, 20]
   ! Stiff   : mu large, e.g. m = 1000, t = [0, 3000]
   !     y'' = m(1 - y^2)y' - y
   ! subject to:
   !     y(0) = 2, y'(0) = 1
   ! Ref: 
   !     https://www.mathworks.com/help/matlab/math/differential-equations.html
   subroutine fvdp(neq,t,u,du)
      implicit none
      integer, intent (in)  ::  neq
      real(8), intent (in)  ::  t, u(neq)
      real(8), intent (out) :: du(neq)
      ! local vars
      real(8) :: y, dy(neq)
      ! unpack y
      y  = u(1)
      dy(1) = u(2)
      ! compute desired d(n)y
      dy(2) = 5d0 * (1 - y*y) * dy(1) - y
      ! assemble du
      du = dy  
   end subroutine fvdp

   ! roots van der Pol
   subroutine frvdp (neq, t, y, ng, gout)
      integer, intent(in)  :: neq, ng
      real(8), intent(in)  :: t, y(3)
      real(8), intent(out) :: gout(2)
      gout(1) = y(1) - 1d0
   end


   !> 1st degree, 3-equations coupled system
   ! Stiff problem from chemical kinetics, consisting of three rate equations:
   !
   !     dy1/dt = -.04*y1 + 1.E4*y2*y3
   !     dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
   !     dy3/dt = 3.E7*y2**2
   !
   ! on the interval t = .4, 4, 4e1, 4e2, .., 4e10
   ! With initial conditions 
   !
   !     y1=1, y2=y3= 0. 
   !
   ! In addition, we want to find the values of t, y1, y2, and y3 at which
   !
   !     (1) y1 reaches the value 1.e-4, and
   !     (2) y3 reaches the value 1.e-2.
   !
   ! Use specific ATOL for each equation (1e6, 1e-10 and 1e-6 respectively). The
   ! jacobian is known.
   !
   subroutine  fchem3 (neq, t, y, dy)
      integer, intent(in)  :: neq
      real(8), intent(in)  :: t, y(3)
      real(8), intent(out) :: dy(3)

      dy(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
      dy(3) = 3.d7*y(2)*y(2)
      dy(2) = -dy(1) - dy(3)      
   end
   !
   ! jacobian
   subroutine  fjchem3 (neq, t, y, ml, mu, pd, nrpd)
      integer  :: neq, ml, mu, nrpd
      real(8)  :: t, y(3), pd(nrpd,3)
      pd(1,1) = -.04d0
      pd(1,2) = 1d4*y(3)
      pd(1,3) = 1d4*y(2)
      pd(2,1) = .04d0
      pd(2,3) = -pd(1,3)
      pd(3,2) = 6d7*y(2)
      pd(2,2) = -pd(1,2) - pd(3,2)
   end
   !
   ! roots (f(u)=0)
   subroutine frchem3 (neq, t, y, ng, gout)
      integer, intent(in)  :: neq, ng
      real(8), intent(in)  :: t, y(3)
      real(8), intent(out) :: gout(2)
      gout(1) = y(1) - 1d-4
      gout(2) = y(3) - 1d-2
   end

end program ode_scifor_test
