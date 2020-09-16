! FFT tests 
! Library: FFTW 3 (http://www.fftw.org)
! FFTW3 - testing
! RSolano, Feb/27/2020
!
! To compile:
! 		gfortran fft_scifor-test.f90 `pkg-config --cflags --libs fftw3` -logpf!

program fft_test

	use fftw3
	use ogpf

	implicit none
	 
	integer, parameter :: N = 1024
	real(8) :: y(N), t(N)
	complex(8) :: yFFT(N/2+1)
	integer(8) :: tPlan
	
	real(8) :: f, p, frange(N/2+1)
	
	type(gpf) :: plt
	integer :: i

	f = 50  ! Hz
	p = 1/f 	! s
	frange = f/2 * [ (i, i=1, N/2+1) ] / (N/2+1)
	
	t = p * [ (i, i=1,N) ]
	y = 1.0*sin(2*3.14*5*t)
	y = y + 2.0*sin(2*3.14*10*t)
	y = y + 0.5*sin(2*3.14*15*t)
	y = y + 0.25*sin(2*3.14*12.5*t)
		
	call plt%plot(t,y,'with lines')
	
	call dfftw_plan_dft_r2c_1d(tPlan, N, y, yFFT, FFTW_ESTIMATE)
	call dfftw_execute_dft_r2c(tPlan, y, yFFT)
	call dfftw_destroy_plan(tPlan)
   
	call plt%plot( frange, abs(yFFT)/N*2,'with lines')
     
end program        