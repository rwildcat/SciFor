! https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/language-reference/a-to-z-reference/o-to-p/pointer-fortran.html

 INTEGER, POINTER :: arrow (:)
 INTEGER, ALLOCATABLE, TARGET :: bullseye (:,:)
 integer, allocatable :: foo(:)

 ! The following statement associates the pointer with an unused
 ! block of memory.

 ALLOCATE (arrow (5), STAT = ierr)
 IF (ierr.eq.0) then
 	WRITE (*,'(/1x,a)') 'ARROW allocated'
 	print *, 'shape =', shape(arrow)
 end if
 	
 arrow = 5
 PRINT *, arrow

! ??
 foo = arrow+3
 print *, 'foo: shape=', shape(foo)
 print *, '=', foo
! ???
 foo = [11,22,arrow, 33,44]
 print *, 'foo: shape=', shape(foo)
 print *, '=', foo



 ALLOCATE (bullseye (8,3), STAT = ierr)
 IF (ierr .eq. 0) then
 	WRITE (*,*) 'BULLSEYE allocated'
 	print *, 'shape =', shape(bullseye)
 end if

 bullseye = 1.
 bullseye (1,:) = 10.
 WRITE (*,'(1x,8i3.0)') bullseye
 print *, bullseye

 ! The following association breaks the association with the first
 ! target, which being unnamed and unassociated with other pointers,
 ! becomes lost. ARROW acquires a new shape.

 arrow => bullseye (2:7,2)
 WRITE (*,'(/1x,a)') 'ARROW is repointed & resized, all the 5s are lost'
 WRITE (*,'(1x,8i8.0)') arrow

 NULLIFY (arrow)
 IF (.NOT.ASSOCIATED(arrow)) WRITE (*,'(/a/)') ' ARROW is not pointed'

 DEALLOCATE (bullseye, STAT = ierr)
 IF (ierr.eq.0) WRITE (*,*) 'Deallocation successful.'

 END