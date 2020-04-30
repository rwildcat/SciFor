! Shortcut for enabling FFTW3
! To compile module:
!	gfortran -c FFTW3.f90 `pkg-config --cflags fftw3`

module FFTW3

    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'

end module