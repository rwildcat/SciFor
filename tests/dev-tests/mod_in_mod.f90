! progrmam to test behavior of a module inside a module
! rsolano, apr/28/2020

module constants
	implicit none
	private
	real, parameter, public :: pi = 3.1415926536  
end module

module areals
	use constants
	implicit none
	private
	public::area_circ, pi
contains	
	real function area_circ(r)	
		real :: r
		area_circ = pi*r*r
	end function
end module


program modtest

	use areals
	implicit none
	real :: cr, carea
	
	cr = 3.0
	carea = area_circ(cr)
	print *, 'pi =', pi
	print *, 'Circle: r=', cr, ', area=', carea
	
end program
	
	
		