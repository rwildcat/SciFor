! SciFor - Interp tests
! Library: FGSL/GSL - GNU Scientific Library
!
! rsolano, feb/27/2020
!
!  $ gfortran interp_scifor-test.f90 `pkg-config --cflags --libs scifor` -l ogpf
!

program interp_scifor
	
	use scifor
	use ogpf
	
	implicit none
	
	type(gpf) :: plt
	integer, parameter :: ndat=6, nplot=100
	real(8) :: xdata(ndat),  ydata(ndat)
	real(8) :: xplot(nplot), yplot(nplot), yplot2(nplot)
	
	! field data
	xdata = [ 1.0, 2.1, 2.9, 3.8, 5.2, 6.4 ]
	ydata = [ 0.5, 3.4, 5.8, 4.1, 1.9, 0.5 ]
	
	! desired interp range
	call numfor_linspace(xplot, xdata(1), xdata(ndat))
	
	! interpolate: linear, poly, cspline, akima, steffen, cspline_p, akima_p
	yplot  = interp_interp(xplot, xdata, ydata, 'poly')
	yplot2 = interp_interp(xplot, xdata, ydata, 'cspline')
	
	! show	
	call plt%plot( &
		xdata, ydata,  'with points ps 3 title "Data"', 'x1y1',			&
		xplot, yplot,  'with lines  lt 1 title "Interp Poly"', 'x1y1', 		&
		xplot, yplot2, 'with lines  lt 2 title "Interp CSplines"', 'x1y1')
		

end program
