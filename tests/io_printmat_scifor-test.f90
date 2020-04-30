! Testing scifor, generic section
! R Solano, Apr/2020
!
! program to test io_printmat_scifor
!
! To compile:
!		# gfortran io_printmat_scifor-test.f90  `pkg-config --cflags --libs scifor` 
!

program test_printmat

   use SciFor
   implicit none

	! make some matrices
   real(8) :: ar8(3,5)
   real(4) :: ar4(3,5)
   integer :: ai(3,5)

	! fill with random contents
   call random_number(ar8)
   call random_number(ar4)
   ai = floor(100*ar4)
   
   ! print double
   call io_print_mat(ar8)
   call io_print_mat(ar8, 'f10.4')
   call io_print_mat(ar8, '3xe10.4')

	! print single
   call io_print_mat(ar4)
   call io_print_mat(ar4, 'f10.4')

	! print integer
   call io_print_mat(ai)
   call io_print_mat(ai, 'i10.8')
   
end program
