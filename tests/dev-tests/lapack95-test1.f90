! https://cyber.dabamos.de/programming/modernfortran/lapack.html
! 
! $ gfortran -L/usr/local/lib -I/usr/local/lib/lapack95_modules \
!   lapack95-test1.f90 -llapack95 -llapack -lblas -ltmglib
!
! $ ./a.out
! Solution (x, y): 1.0000, 3.0000 
! 

! lapack95
program main

    use :: la_precision, only: wp => dp
    use :: f95_lapack,   only: la_gesv
    
    implicit none
    real(kind=wp) :: a(2, 2)
    real(kind=wp) :: b(2)

    a(2, 2) = reshape([ 2., 3., 1., 1. ], [ 2, 2 ])
    b(2)    = [ 5., 6. ]

    call la_gesv(a, b)

    print '(a, 2(f0.4, ", "))', 'Solution (x, y): ', b

end program main