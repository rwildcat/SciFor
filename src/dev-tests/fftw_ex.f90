! example from: appleparan/fftw_ex.f90
! 		https://gist.github.com/appleparan/c048c44668ede7ef28ba63c660b6dcf3

program main

	use ISO_C_BINDING
	implicit none

	include 'fftw3.f03'

	real(8), allocatable :: xin(:)
	complex(8), allocatable :: xout(:)
	
	integer(8) :: plan_c2r, plan_r2c
	integer :: n, i

	real(8), allocatable :: r(:), corr1(:), corr2(:)
	integer(8) :: plan_c2r_corr(2), plan_r2c_corr(2)
	integer :: n_corr, i_corr

	n = 9
	allocate( xin(n) )
	allocate( xout(n) )

	! dummy content for planning
	xin = 0
	xout = 0

	! sound like a plan
	call dfftw_plan_dft_r2c_1d(plan_r2c, n, xin, xout, FFTW_ESTIMATE)
	!call dfftw_plan_dft_c2r_1d(plan_c2r, n, xin, xout, FFTW_ESTIMATE)

	! fill real data to run the plan
! 	xin = [1., 2., 3., 4., 5.]
	xin = [1,0,0,0,0,0,0,0,0]
!  fft() =  {0.333333 x  9 times}

	! fft in-place
	print *, "IN: "
	print *, xin

	call dfftw_execute_dft_r2c(plan_r2c, xin, xout)
	
	print *, "IN: "
	print *, xin

	print *, "OUT: "	
	print *, abs(xout)

	stop
	
	
	call dfftw_execute_dft_c2r(plan_c2r, xin, xout)

	do i = 1, n
	  ! normalization
	  xin(i) = xin(i) / n
	end do

	print *, "IN2: "
	write(*, *) (xin(i),i = 1,n)

	call dfftw_destroy_plan(plan_r2c)
	call dfftw_destroy_plan(plan_c2r)

	deallocate(xin)
	deallocate(xout)


	! corr
	n_corr = 192
	allocate(r(n_corr))
	allocate(corr1(n_corr+2))
	allocate(corr2(n_corr+2))

	call dfftw_plan_dft_r2c_1d(plan_r2c_corr(1), n_corr, corr1, corr1, FFTW_MEASURE)
	call dfftw_plan_dft_c2r_1d(plan_c2r_corr(1), n_corr, corr1, corr1, FFTW_MEASURE)
	call dfftw_plan_dft_r2c_1d(plan_r2c_corr(2), n_corr, corr2, corr2, FFTW_MEASURE)
	call dfftw_plan_dft_c2r_1d(plan_c2r_corr(2), n_corr, corr2, corr2, FFTW_MEASURE)

	open(100, FILE='corr.plt', ACTION='READ')

	do i = 1, n_corr
	  read(100, *) r(i),corr1(i), corr2(i)
	  print *, r(i), corr1(i), corr2(i)
	end do
	close(100)
	call dfftw_execute_dft_r2c(plan_r2c_corr(1), corr1, corr1)
	call dfftw_execute_dft_r2c(plan_r2c_corr(2), corr2, corr2)

	open(100, FILE='E1D.plt', ACTION='WRITE')
	do i = 1, n_corr/2
	  write(100, '(I4, 3ES20.9)') &
			i,  &
			2.0*corr1(2*(i-1)+1), &
			2.0*corr2(2*(i-1)+1)

	end do
	close(100)

	deallocate(r)
	deallocate(corr1)
	deallocate(corr2)
	do i = 1, 2
	  call dfftw_destroy_plan(plan_r2c_corr(i))
	  call dfftw_destroy_plan(plan_c2r_corr(i))
	end do
	 
end program