! https://cyber.dabamos.de/programming/modernfortran/lapack.html
!
! To compile:
! 		$ gfortran -L/usr/local/lib/ lapack-test1.f90 -llapack -lblas
!		$ gfortran lapack-test1.f90 -framework accelerate
!
! 		$ ./a.out
! Solution (x, y): 1.0000, 3.0000

! lapack
program main

    implicit none
    external :: sgesv
    real    :: a(2, 2)  ! Matrix A.
    real    :: b(2)     ! Matrix B.
    real    :: pivot(2) ! Pivot indices (list of swap operations).
    integer :: rc       ! Return code.

    a = reshape([ 2., 3., 1., 1. ], [ 2, 2 ])
    b = [ 5., 6. ]

    call sgesv(2, 1, a, 2, pivot, b, 2, rc)

    if (rc == 0) then
        print '(a, 2(f0.4, ", "))', 'Solution (x, y): ', b
    else
        print '(a, i0)', 'Error: ', rc
    end if
    
end program main