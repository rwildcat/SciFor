!> Program to test `scifor_ode`
!
! To compile:
!  `$ gfortran ode_scifor-test.f90 `pkg-config --cflags --libs scifor`

!
program ode_scifor_test

   !use scifor_ode
   use scifor
   implicit none
   
   type(Odesolv) :: ode 

   ! example 1: 
   !     1st deg, stiff + roots: flame radius, y' = y^2 - y^3
   !
	call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3] )
   ! call ode%solve( fdy, [2d-3], trng=[0d0, 2d0/2d-3], froots=fgex, nfroots=1 )
   ! call ode%solve( fdy, [2d-3], tpts=[0d0, 1d1, 3d2, 4.5d2, 4.75d2, 4.9d2, 5d2, 5.1d2, 6d2, 10d2] )
   ! call ode%solve( fdy, [2d-3], tpts=[0d0, 1d1, 3d2, 4.5d2, 4.75d2, 4.9d2, 5d2, 5.1d2, 6d2, 10d2], froots=fgex, nfroots=1 )

   ! example 2: 
   !     2nd deg, stiff/no stiff: van der Pool oscillator, y'' = m(1 - y^2)y' - y
   !
   ! call ode%solve( fvdp, [2d0,0d0], trng=[0d0, 3000d0] )

   ! example 3: 
   !     1st deg, 3-eq system + roots + jacobian: 3 substances ratios
   !
   ! call ode%solve(                &
   !    fchem3, [1d0, 0d0, 0d0],    &
   !    tpts=[0d0, 4d-1, 4d0, 4d1, 4d2, 4d3, 4d4, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10],  &
   !    rtol=[1d-4], atol=[1.d-6, 1.d-10, 1.d-6],   &
   !    jac=fjchem3, froots=frchem3, nfroots=2 )
   
contains

   !> 1st degree, stiff ODE system - Flame radius 
   !		y' = y^2 - y^3
   ! subject to: 
   !		y(0) = delta
   !		0 <= t <= 2/delta
   subroutine fdy(neq,t,u,du)
      implicit none
      integer, intent (in)  ::  neq
      real(8), intent (in)  ::  t, u(neq)
      real(8), intent (out) :: du(neq)
      ! local vars
      real(8) :: y, dy(neq)
      ! unpack y
      y  = u(1)
      ! compute desired d(n)y
      dy(1) = y**2 - y**3
      ! assemble du
      du = dy  
   end subroutine fdy
   !
   subroutine fgex (neq, t, y, ng, gout)
      integer :: neq, ng
      real(8) :: t, y(1), gout(1)
      gout(1) = y(1) - 5d-2
   end subroutine


   !> 2nd degree, no/yes stiff ODE system - van der Pool oscillator
   !
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
      dy(2) = 1000d0 * (1 - y*y) * dy(1) - y
      ! assemble du
      du = dy  
   end subroutine fvdp


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
   subroutine  fchem3 (neq, t, y, ydot)
      integer  :: neq
      real(8)  :: t, y(3), ydot(3)
      ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
      ydot(3) = 3.d7*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)      
   end
!
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
   subroutine frchem3 (neq, t, y, ng, gout)
      real(8)  :: t, y(3), gout(2)
      integer  :: neq, ng
      gout(1) = y(1) - 1d-4
      gout(2) = y(3) - 1d-2
   end



end program ode_scifor_test
