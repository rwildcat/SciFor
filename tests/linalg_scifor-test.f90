! Testing SciFor (LAPACK)
! rsolano, feb/25/2020
!
! Test may rise:
! 	"Note: The following floating-point exceptions are signalling: IEEE_DENORMAL"
!
! To remove warning:
! 	- Use cflag fpe-summary=invalid,zero,overflow
!
! For example:
!	0. $ export FERR=fpe-summary=invalid,zero,overflow
! 	1. $ gfortran SciFor.f90 -f$FFERR scifor-test.f90 -framework Accelerate
! 	2. $ gfortran SciFor.f90 -f$FFERR scifor-test.f90 -l lapack
!
! To compile:
!		$ gfortran linalg_scifor-test.f90 `pkg-config --cflags --libs scifor` 
!

program linalg1

   use scifor
   implicit none
   
   integer, parameter :: nmax=3, nRHS=2
   real(8) :: A(nmax,nmax)
   real(8) :: b(nmax)      , x(nmax)
   real(8) :: bn(nmax,nRHS), xn(nmax,nRHS)
   integer :: n
   
   ! system size
   n = nmax

   print *, 'Testing nRHS=2'

   ! set data   
   A = reshape( [ 3,2,0 , 1,-1,0 , 0,5,1], [3,3], order=[2,1] )
   bn = reshape( [2,4,-1, 2,4,1], [3,2] )
   b = [2,4,1]   

   ! solve A*bn = xn
   xn = linalg_solve(A, bn)
               
   ! show 
   call io_print_mat(A)
   print *, ''
   call io_print_mat(xn)
   
   ! ---
   
    print *,''
    print *, 'Test nRHS=1'
    
	! solve A*b = x
   x = linalg_solve(A, b)
               
   ! show 
   call io_print_mat(A)
   print *, ''
   print *, x

end program randomsys1
