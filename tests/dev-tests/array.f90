program array

	implicit none
	
	integer :: i, j, a(12)
	
	a = [(i, i=1,size(a))]
	
	call print_matd(1d0*reshape(a, [2,6]), 'ES10.2')

contains 

subroutine print_matd(A, fmt_)
      real(8), intent(in) :: A(:,:)
      character(*),optional, intent(in) :: fmt_
      
      ! my vars
      character(30) :: fmt = 'ES12.4'
      integer :: i
      
      ! fmt_ = '*' --> ''
      if (present(fmt_)) then
         if (trim(fmt_) == '*') then            
            fmt = ''
         else
            fmt = fmt_
         end if         
      end if
      
      if (len_trim(fmt) > 0) then
         fmt = '(99(' // trim(fmt) // '))' 
      end if
            
      do i=1,size(A,1)
         if (len_trim(fmt) > 0 ) then
            print fmt, A(i,:)
         else
            print *, A(i,:)
         end if
      end do
      
   end subroutine
    


end program
