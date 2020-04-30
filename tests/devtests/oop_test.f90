module namespacename ! the module name defines the namespace

		type classname ! classname is the class prototype name
			integer :: stat_field1
			complex         :: inst_field1
			
		contains
			! declare a static constructor        
			procedure classname_cctor
			
			! declare some instance constructors
			procedure classname_ctor0
			procedure classname_ctor1
			procedure classname_ctor2
		end type classname   

	contains

		! implement static constructor        
		subroutine classname_cctor(this)
			class(classname) :: this
			integer :: i
			stat_field1 = -1
		end subroutine

		! implement instance constructor with no arguments
		subroutine classname_ctor0(this)
		  type(classname) :: this
		  this%inst_field1 = cmplx(0.,0.) 
		end subroutine classname_ctor0

		! implement instance constructor with one argument
		subroutine classname_ctor1(this,Value)
		  type(classname) :: this
		  real,value :: Value
		  this%inst_field1 = cmplx(Value,0.)
		end subroutine classname_ctor1

		! implement instance constructor with two arguments
		subroutine classname_ctor2(this,Value1,Value2)
		  type(classname) :: this
		  real,value :: Value1,Value2
		  this%inst_field1 = cmplx(Value1,Value2)
		end subroutine classname_ctor2

end module namespacename
	