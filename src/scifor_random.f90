!> # Unified interface to `randlib` (randlib90) library 
! Ramon Solano
! ramon.solano @ gmail.com
! COlima, Mexico
!
!
! ## RANDLIB - Library of Routines for Random Number Generation
! Routines are provided that return random deviates from the following distributions:
!
! * Real uniform
! * Integer uniform
! * Beta
! * Chi-square
! * Exponential
! * Gamma
! * Multivariate normal
! * Non-central chi-squared
! * Non-central F
! * Univariate normal
! * Random permutations
! * Binomial
! * Negative Binomial
! * Multinomial
! * Poisson

!
! Barry W. Brown
! Departments of Biomathematics and Biostatistics 
! Software developers: John Venier, Dan Serachitopol
! * [randlib90](https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/27)
!
module scifor_random_mod

   use user_set_generator
   use random_standard_uniform_mod  ! nfdist 1
   use random_normal_mod            ! nfdist 2
   use random_beta_mod              ! nfdist 3
   use random_binomial_mod          ! nfdist 4
   use random_chisq_mod             ! nfdist 5
   use random_exponential_mod       ! nfdist 6

   implicit none
   
   private

   type, public :: Scifor_random
      private
      character(20)  :: fdists(10)   = [ character(20) ::   &  ! fdists names
         'uniform','normal','beta','binomial','chisq','exp']      
      character(20)  :: fdist        = 'uniform'    ! default
      integer        :: nfdist       = 1            ! fdist number
      integer        :: seed         = 0            ! init seed
      real           :: params(3)    = 0            ! fdist params
      integer        :: nparams      = 0            ! n of params
      logical        :: initialized  = .False.      ! already init?

   contains
      private      
      procedure            :: rand_rd0    ! 0-d real
      procedure            :: rand_rd1    ! 1-d real
      procedure            :: rand_rd2    ! 2-d real

      procedure, public    :: init
      procedure, public    :: set_seed
      procedure, public    :: info
      generic, public      :: rand => rand_rd0, rand_rd1, rand_rd2


   end type

   ! optional constructor
   interface Scifor_random
      module procedure constructor
   end interface

contains

!> Constructor (opt)
! Provides a canonical constructor function to be called as e.g. `ran = Scifor_random()`
! However, it is not called intrinsecally when the object is constructed, but it is up to 
! the user to call it.
! To work around, the function `*%init()` is provided and is called by the
! constructor as only task. This way, `init()` can be called when necessary.
! Constructor kept only for style completness
!
function constructor(fdist, params, seed) result(Ran)
   type(Scifor_random) :: Ran
   character(*), optional, intent(in) :: fdist    ! function distribution
   real, optional, intent(in)         :: params(:)   ! function parameters
   integer, optional,intent(in)       :: seed        ! random seed

   ! print *, '>> Scifor_random%constructor() called'

   if (present(fdist)) Ran%fdist = fdist
   if (present(seed))   Ran%seed = seed

   if (present(params)) then
      Ran%params = params
      Ran%nparams = size(params)
   end if

   ! init
   call Ran%init(fdist, params, seed)

end function constructor

! ---

! initialization
subroutine init(this, fdist, params, seed)
   class(Scifor_random)                :: this
   character(*), optional, intent(in)  :: fdist       ! function distribution name
   real, optional, intent(in)          :: params(:)   ! function parameters
   integer, optional, intent(in)       :: seed        ! random seed
   ! my vars
   character                           :: numstr(80)
   
   ! print *, '>> Scifor_random%init() called'

   ! defaults already set

   ! set new optional defaults
   if (present(fdist))  this%fdist  = fdist
   if (present(seed))   this%seed   =  seed

   if (present(params)) then
      this%params = params
      this%nparams = size(params)
   end if

   ! set nfdist
   select case (this%fdist)

      case ('uniform')
         this%nfdist = 1
         if (this%nparams > 0) then
            error stop '** Error: *uniform* distribution: please specify no params'
         end if

      case ('normal')
         this%nfdist = 2
         ! params: (m, sd); m=mean, sd=std dev
         if (this%nparams == 0) then
            this%nparams = 2
            this%params(1:2) = [0,1]
         else if (this%nparams /= 2) then
            error stop '** Error: *normal* distribution: please specify either 0 or 2 params'
         end if 

      case ('beta')
         this%nfdist = 3
         ! params: (a, b); a=b, b=b
         if (this%nparams == 0) then
            this%nparams = 2
            this%params(1:2) = [5.0 ,5.0]
         else if (this%nparams /= 2) then
            error stop '** Error: *beta* distribution: please specify either 0 or 2 params' 
         end if 

      case ('binomial')
         this%nfdist = 4
         ! params: (n, p); n=n of events, p=prob of each event
         if (this%nparams == 0) then
            this%nparams = 2
            this%params(1:2) = [1.0, 0.5]
         else if (this%nparams /= 2) then
            error stop '** Error: *binomial* distribution: please specify either 0 or 2 params'
         end if 
         
      case ('chisq')
         this%nfdist = 5
         ! params: (df); df=degrees of freedom
         if (this%nparams == 0) then
            this%nparams = 1
            this%params(1) = 1.0
         else if (this%nparams /= 1) then
            error stop '** Error: *chisq* distribution: please specify either 0 or 1 param'
         end if 

      case ('exp')
         this%nfdist = 6
         ! params: (m); m=mean; beta=1/lambda
         if (this%nparams == 0) then
            this%nparams = 1
            this%params(1) = 1.0
         else if (this%nparams /= 1) then
            error stop '** Error: *exp* distribution: please specify either 0 or 1 param'
         end if 

      case default
         error stop '** Error: unknown random function distribution : "' // trim(this%fdist) // '"'
   end select

   ! intialized as done 
   this%initialized = .True.

   ! set seed
   call this%set_seed(this%seed)
   
end subroutine

! ---

subroutine set_seed(this, seed)
   class(Scifor_random)          :: this
   integer, optional, intent(in) :: seed        ! random seed
   ! vars
   character(80)                 :: phrase      ! phrase to init seed

    ! init if required
   if (.not. this%initialized) call this%init()

   if (present(seed)) then
      this%seed = seed
   end if

   if (this%seed == 0) this%seed = 1

   ! init accordingly
   if (this%seed == 1) then
      ! print *, '>> seeding: seed=1'
      call time_set_seeds()
   else
      write(phrase,*) this%seed
      ! print *, '>> seeding: phrase = [' // trim(phrase) // ']'
      call phrase_set_seeds(phrase)
   end if

   
end subroutine

! ---

subroutine info(this)
   class(Scifor_random) :: this
   print *, '>> fdist=' // trim(this%fdist) // ', seed=', this%seed, ', params=', this%params(:this%nparams)
end subroutine

! ---

! 0-d array
subroutine rand_rd0(this, a)
   class(Scifor_random)               :: this
   real, intent(out)                  :: a
  
   ! print *,  '>> rand_rd0() called'

   ! init if required
   if (.not. this%initialized) call this%init()

   ! generate accordingly
   select case (this%nfdist)
      case (1)
         a = random_standard_uniform()
      case (2)
         a = random_normal(this%params(1), this%params(2))
      case (3)
         a = random_beta(this%params(1), this%params(2))
      case(4)
         a = random_binomial(int(this%params(1)), this%params(2))
      case(5)
         a = random_chisq(this%params(1))
      case(6)
         a = random_exponential(this%params(1))
      case default   ! error
         error stop '** Error: random function distribution unknown'
   end select

end subroutine

! ---

! 1-d array
subroutine rand_rd1(this, a)
   class(Scifor_random)               :: this
   real, intent(out)                  :: a(:)
   ! local vars
   integer              :: i

   ! print *,  '>> rand_rd1() called'

   ! init if required
   if (.not. this%initialized) call this%init()

   ! generate accordingly
   select case (this%nfdist)
      case (1)
         a = [ (random_standard_uniform(), i=1, size(a)) ]
      case (2)
         a = [ (random_normal(this%params(1), this%params(2)), i=1, size(a)) ]
      case (3)
         a = [ (random_beta(this%params(1), this%params(2)), i=1, size(a)) ]
      case(4)
         a = [ (random_binomial(int(this%params(1)), this%params(2)), i=1, size(a)) ]
      case(5)
         a = [ (random_chisq(this%params(1)), i=1, size(a)) ]
      case(6)
         a = [ (random_exponential(this%params(1)), i=1, size(a)) ]
      case default   ! error
         error stop '** Error: random function distribution unknown'
   end select

end subroutine

! ---

! 2-d array
subroutine rand_rd2(this, a)
   class(Scifor_random) :: this
   real                 :: a(:,:)
   ! local vars
   integer              :: i

   ! print *,  '>> rand_rd1() called'

   ! init if required
   if (.not. this%initialized) call this%init()

   ! generate accordingly
   select case (this%nfdist)
      case (1)
         a = reshape([ (random_standard_uniform(), i=1, size(a)) ], shape=shape(a))
      case (2)
         a = reshape([ (random_normal(this%params(1), this%params(2)), i=1, size(a)) ], shape=shape(a))
      case (3)
         a = reshape([ (random_beta(this%params(1), this%params(2)), i=1, size(a)) ], shape=shape(a))
      case(4)
         a = reshape([ (random_binomial(int(this%params(1)), this%params(2)), i=1, size(a)) ], shape=shape(a))
      case(5)
         a = reshape([ (random_chisq(this%params(1)), i=1, size(a)) ], shape=shape(a))
      case(6)
         a = reshape([ (random_exponential(this%params(1)), i=1, size(a)) ], shape=shape(a))
      case default   ! error
         error stop '** Error: random function distribution unknown'
   end select


   print *,  '>> rand_rd2() called'
end subroutine



end module
