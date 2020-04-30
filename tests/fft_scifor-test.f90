! SciFor tests
! Library: FFTW 3 (http://www.fftw.org)
!
! RSolano, Feb/27/2020
!
! To compile the SciFor:
!		$ gfortran fft_scifor-test.f90 -l ogpf `pkg-config --cflags --libs scifor` 
!

program fft_scifor_test

	use scifor
	use ogpf

	implicit none
	 
	integer, parameter :: N = 1024 * 1
	real(8) :: y(N), t(N), yi(N), rnoise(N)
	complex(8) :: yFFT(N)
	
	real(8) :: f, p, frange(N/2+1)
	
	type(gpf) :: plt
	integer :: i

	f = 100  ! Hz
	p = 1/f 	! s
	
	! frange <0 - f/2] Hz, Nyquist!
	frange = f/2 * [ (i, i=1, N/2+1) ] / (N/2+1)
	
	! make a composite signal
	t = p * [ (i, i=1,N) ]
	y =     1.5*sin(2*3.14 * 5* t)
	y = y + 1.0*sin(2*3.14 *10* t)
	y = y + 0.5*sin(2*3.14 *20* t)

	! add some noise
	call random_number(rnoise)
	y = y + 4*rnoise - 2

	! fft
	yFFT = fft_fft(y)  
	
	! ifft
	yi = fft_ifft(yFFT)

	! max
	print *, 'err =', maxval(abs(yi-y))

	! now a nice plot
	call plt%multiplot(3,1) 

	call plt%title('Original data (y = 1.5 sin(5) + sin(10t)+ .5 sin(20t))')
	call plt%xlabel('t (s)')
	call plt%ylabel('A')
	call plt%plot(t, y, 'w lines')

	call plt%title('Transformed data (FFT)')
	call plt%xlabel('f (Hz)')
	call plt%ylabel('PS')
	call plt%plot(frange, abs(yFFT(1:n/2+1)), 'w lines')
	
	call plt%title('Reconstructed data (IFFT)')
	call plt%xlabel('t (s)')
	call plt%ylabel('A')
	call plt%plot(t, yi, 'w lines')

	
	  
end program        