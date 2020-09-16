! testing a dynamically growing array
! rsolano, aug/11/2020

program dyngrow

   implicit none
   real, allocatable, target :: a1(:), a2(:), a(:)
   real, pointer        :: p(:)
   integer              :: i, nchunk
   logical              :: doAnother
   real                 :: r

   nchunk = 5
   allocate(a1(nchunk))
   p => a1

   doAnother = .True.
   i = 0

   do while (doAnother)

      call random_number(r)
      print *, 'i=', i+1, ', r=', r

      if (r > 0.5) then
         print *, '--'
         i = i + 1

         ! grow array if needed
         if (i>size(p)) then
            print *, '** p full, i=', i,' size=', size(p)
            if (associated(p,target=a1)) then
               print *, '  - a1 full; allocating into a2'
               allocate( a2(size(a1)+nchunk) )  ! make new room in a2
               a2(:size(a1)) = a1               ! copy a1 into a2
               p => a2                          ! point p to a2
               deallocate(a1)                   ! release a1
            else
               print *, '  - a2 full; allocating into a1'
               allocate( a1(size(a2)+nchunk) )  ! make new room in a1
               a1(:size(a2)) = a2               ! copy a1 into a2
               p => a1                          ! point p to a2
               deallocate(a2)                   ! release a1
            end if
         end if

         ! assign
         p(i) = i+r
         doAnother = r < 0.95
      end if

   end do

   ! set result
   a = p(:i)

   ! error? -- yes!
   !deallocate(p)

   ! do it right
   if (allocated(a1)) deallocate(a1)
   if (allocated(a2)) deallocate(a2)

   print '(/,A)', 'end loop'

   print *, 'i =', i
   print *, 'size(p), p()    =', size(p), p
   print *, 'allocated(a1) :', allocated(a1)
   print *, 'allocated(a2) :', allocated(a2)
   print *, ''
   print *, 'size(a), a()    =', size(a), a

   
end program dyngrow