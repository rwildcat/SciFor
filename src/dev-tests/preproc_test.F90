! program to test conditional compiling w/ fortran preprocessor
! to compile:
!		$ gfortran -cpp  -D UOS preproc_test.F90 
!

program test_Dopt

	implicit none
	character (len=1) :: pathsep
	
	pathsep = "?"
	
#	ifdef WOS
		pathsep = "\"
#	endif

#	ifdef UOS
		pathsep = "/"
#	endif

	write (*, '( "pathsep is >", A1, "<")' )  pathsep

end program test_Dopt