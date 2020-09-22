! Testing scifor, random number section
! R Solano, Sep/2020
!
! program to test scifor_random()
!
! To compile:
!		$ gfortran random_scifor-test.f90  `pkg-config --cflags --libs scifor` 
!

program test_random

   use scifor
   implicit none

   ! my vars
   type(Scifor_random)  :: rgen
   integer, parameter   :: n=3
   real(4)              :: x, xv(n*n), yv(n), xm(n,n)
   integer              :: i

   ! testing constructor
   ! rgen = Scifor_random(fdist='normal', params=[10.0, 1.0], seed=0)

   ! testing init 
   ! call rgen%init()
   ! call rgen%init(fdist='normal', params=[10.0, 1.0], seed=0)

   ! init w/ default and custom params
   ! call rgen%init('normal') 
   ! call rgen%init('normal', [10.0, 3.0])

   ! init w/ default and custom params
   ! call rgen%init('beta')
   ! call rgen%init('beta', [8.0, 2.0])

   ! init w/ default and custom params
   ! call rgen%init('binomial')
   ! call rgen%init('binomial', [10.0, 0.5])

   ! init w/ default and custom params
   ! call rgen%init('chisq')
   ! call rgen%init('chisq', [3.0])

   ! exponential
   ! call rgen%init('exp')
   ! call rgen%init('exp', [5.0])

   ! some data, using real 0-d argument
   ! call rgen%set_seed(0)

   call rgen%init('beta', [8.0,2.0], 2)
   
   do i=1, n*n
      call rgen%rand(x)
      print '(F10.5)', x
   end do
   call rgen%info

   call rgen%init(seed=2)
   call rgen%rand(xv)
   print *, xv
   call rgen%info

   call rgen%init(seed=2)
   call rgen%rand(xm)
   call io_print_mat(xm)
   call rgen%info


 
end program
