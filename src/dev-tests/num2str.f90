! program to test converting nums to strings
program num2str
   implicit none
   character(80) :: fmt,fmt0,nc
   real          :: a(3,4)
   integer        :: i

   fmt =  'es12.4'

   call random_number(a)

   ! assemble fmt
   write(nc,*) size(a,2)
   nc = adjustl(nc)

   write (fmt0, *)  '(' // trim(nc) // '(' // trim(fmt) // '))'
   print *, 'fmt =', fmt0

   print *, 'a ='

   do i=1, size(a,1)
      print fmt0, a(i,:)
   end do


end program
