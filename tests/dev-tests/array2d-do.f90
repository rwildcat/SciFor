program array2d

	implicit none
	integer, parameter :: nr=2, nc=3
	real	:: a(nr,nc)
	integer :: i,j
	
	a(:,1) = [1,2]
! 	a = reshape([1,2, 3,4, 5,6], [nr,nc])
	a = reshape([(i, i=1,nr*nc)], [nr,nc])
	
	print *, 'rank(a) =', rank(a)
	print *, 'shape(a)=', shape(a)
	print *, 'size(a) =', size(a)
	
	print *, 'a='
	do i=1, nr
		print *, a(i,:)
	end do
	
	
end program
	