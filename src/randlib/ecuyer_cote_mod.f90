    MODULE ecuyer_cote_mod
! .. Implicit None Statement ..
      IMPLICIT NONE
! ..
! .. Default Accessibility ..
      PRIVATE
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: a1 = 40014, a1vw = 2082007225, &
        a1w = 1033780774, a2 = 40692, a2vw = 784306273, a2w = 1494757890, &
        default_seed1 = 1234567890, default_seed2 = 123456789, &
        m1 = 2147483563, m2 = 2147483399, n_generators = 32
! ..
! .. Local Scalars ..
      INTEGER :: current_generator = 1
      LOGICAL :: generator_set = .FALSE., seeds_set = .FALSE.
! ..
! .. Local Arrays ..
      INTEGER :: cg1(n_generators), cg2(n_generators), ig1(n_generators), &
        ig2(n_generators), lg1(n_generators), lg2(n_generators)
      LOGICAL :: qanti(n_generators) = .FALSE.
! ..
! .. Public Statements ..
      PUBLIC :: advance_state, get_current_generator, get_seeds, &
        random_large_integer, reinitialize_current_generator, &
        set_all_seeds, set_antithetic, set_current_generator, &
        set_current_seed
! ..
    CONTAINS

!*********************************************************************

      SUBROUTINE advance_state(k)
!----------------------------------------------------------------------
!     SUBROUTINE ADVNST(K)
!               ADV-a-N-ce ST-ate
!     Advances the state  of  the current  generator  by 2^K values  and
!     resets the initial seed to that value.
!     This is  a  transcription from   Pascal to  Fortran    of  routine
!     Advance_State from the paper
!     L'Ecuyer, P. and  Cote, S. "Implementing  a  Random Number Package
!     with  Splitting   Facilities."  ACM  Transactions  on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     K -> The generator is advanced by2^K values
!                                   INTEGER K
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: k
! ..
! .. Local Scalars ..
        INTEGER :: g, i, ib1, ib2
! ..
! .. Executable Statements ..

        IF ( .NOT. seeds_set) CALL set_default

        CALL get_current_generator(g)

        ib1 = a1
        ib2 = a2
        DO i = 1, k
          ib1 = multiply_modulo(ib1,ib1,m1)
          ib2 = multiply_modulo(ib2,ib2,m2)
        END DO
        CALL set_current_seed(multiply_modulo(ib1,cg1(g),m1), &
          multiply_modulo(ib2,cg2(g),m2))

!     NOW, IB1 = A1**K AND IB2 = A2**K

        RETURN

      END SUBROUTINE advance_state

!*********************************************************************

      SUBROUTINE get_current_generator(g)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (OUT) :: g
! ..
! .. Executable Statements ..

        IF ( .NOT. generator_set) THEN
          current_generator = 1
          generator_set = .TRUE.
        END IF

        g = current_generator


        RETURN

      END SUBROUTINE get_current_generator

!*********************************************************************

      SUBROUTINE get_seeds(iseed1,iseed2)
!----------------------------------------------------------------------
!     SUBROUTINE GETSD(G,ISEED1,ISEED2)
!               GET SeeD
!     Returns the value of two integer seeds of the current generator
!     This  is   a  transcription from  Pascal   to  Fortran  of routine
!     Get_State from the paper
!     L'Ecuyer, P. and  Cote,  S. "Implementing a Random Number  Package
!     with   Splitting Facilities."  ACM  Transactions   on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     ISEED1 <- First integer seed of generator G
!                                   INTEGER ISEED1
!     ISEED2 <- Second integer seed of generator G
!                                   INTEGER ISEED1
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (OUT) :: iseed1, iseed2
! ..
! .. Local Scalars ..
        INTEGER :: g
! ..
! .. Executable Statements ..

        IF ( .NOT. seeds_set) CALL set_default

        CALL get_current_generator(g)
        iseed1 = cg1(g)
        iseed2 = cg2(g)
        RETURN

      END SUBROUTINE get_seeds

!*********************************************************************

      FUNCTION multiply_modulo(a,s,m)
!----------------------------------------------------------------------
!     INTEGER FUNCTION MULTIPLY_MODULO(A,S,M)
!                    Returns (A*S) MOD M
!     This is a transcription from Pascal to Fortran of routine
!     MULtMod_Decompos from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     A, S, M  -->
!                         INTEGER A,S,M
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Function Return Value ..
        INTEGER :: multiply_modulo
! ..
! .. Parameters ..
        INTEGER, PARAMETER :: h = 32768
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: a, m, s
! ..
! .. Local Scalars ..
        INTEGER :: a0, a1, k, p, q, qh, rh
! ..
! .. Executable Statements ..

!     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
!      machine. On a different machine recompute H
        IF (a<=0 .OR. a>=m .OR. s<=0 .OR. s>=m) THEN
          WRITE (*,*) ' A, M, S out of order in MLTMOD - ABORT!'
          WRITE (*,*) ' A = ', a, ' S = ', s, ' M = ', m
          WRITE (*,*) ' MLTMOD requires: 0 < A < M; 0 < S < M'
          STOP ' A, M, S out of order in MLTMOD - ABORT!'
        END IF

        IF (a<h) THEN
          a0 = a
          p = 0
        ELSE

          a1 = a/h
          a0 = a - h*a1
          qh = m/h
          rh = m - h*qh
          IF (a1>=h) THEN
            a1 = a1 - h
            k = s/qh
            p = h*(s-k*qh) - k*rh
            DO WHILE (p<0)
              p = p + m
            END DO

          ELSE

            p = 0
          END IF

!     P = (A2*S*H)MOD M

          IF (a1/=0) THEN
            q = m/a1
            k = s/q
            p = p - k*(m-a1*q)
            IF (p>0) p = p - m
            p = p + a1*(s-k*q)
            DO WHILE (p<0)
              p = p + m
            END DO

          END IF
          k = p/qh

!     P = ((A2*H + A1)*S)MOD M

          p = h*(p-k*qh) - k*rh
          DO WHILE (p<0)
            p = p + m
          END DO

        END IF
        IF (a0/=0) THEN

!     P = ((A2*H + A1)*H*S)MOD M

          q = m/a0
          k = s/q
          p = p - k*(m-a0*q)
          IF (p>0) p = p - m
          p = p + a0*(s-k*q)
          DO WHILE (p<0)
            p = p + m
          END DO

        END IF
        multiply_modulo = p

        RETURN

      END FUNCTION multiply_modulo

!*********************************************************************

      FUNCTION random_large_integer()
!----------------------------------------------------------------------
!     INTEGER FUNCTION random_large_integer()
!     Returns a random integer following a uniform distribution over
!     (1, 2147483562) using the current generator.
!     This is a transcription from Pascal to Fortran of routine
!     Random from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Function Return Value ..
        INTEGER :: random_large_integer
! ..
! .. Local Scalars ..
        INTEGER :: curntg, k, s1, s2, z
! ..
! .. Executable Statements ..

!     Get Current Generator
        IF ( .NOT. seeds_set) CALL set_default
        CALL get_current_generator(curntg)
        s1 = cg1(curntg)
        s2 = cg2(curntg)
        k = s1/53668
        s1 = a1*(s1-k*53668) - k*12211
        IF (s1<0) s1 = s1 + m1
        k = s2/52774
        s2 = a2*(s2-k*52774) - k*3791
        IF (s2<0) s2 = s2 + m2
        cg1(curntg) = s1
        cg2(curntg) = s2
        z = s1 - s2
        IF (z<1) z = z + m1 - 1
        IF (qanti(curntg)) z = m1 - z
        random_large_integer = z
        RETURN

      END FUNCTION random_large_integer

!*********************************************************************

      SUBROUTINE reinitialize_current_generator(isdtyp)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: isdtyp
! ..
! .. Local Scalars ..
        INTEGER :: g
! ..
! .. Executable Statements ..

!----------------------------------------------------------------------
!     SUBROUTINE REINITIALIZE_CURRENT_GENERATOR(ISDTYP)
!          INIT-ialize current G-e-N-erator
!     Reinitializes the state of the current generator
!     This is a transcription from Pascal to Fortran of routine
!     Init_Generator from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     ISDTYP -> The state to which the generator is to be set
!          ISDTYP = -1  => sets the seeds to their initial value
!          ISDTYP =  0  => sets the seeds to the first value of
!                          the current block
!          ISDTYP =  1  => sets the seeds to the first value of
!                          the next block
!                                   INTEGER ISDTYP
!----------------------------------------------------------------------
        IF ( .NOT. seeds_set) CALL set_default
        CALL get_current_generator(g)
        SELECT CASE (isdtyp)
        CASE (-1)
          lg1(g) = ig1(g)
          lg2(g) = ig2(g)
        CASE (0)
          CONTINUE
        CASE (1)
          lg1(g) = multiply_modulo(a1w,lg1(g),m1)
          lg2(g) = multiply_modulo(a2w,lg2(g),m2)
        CASE DEFAULT
          STOP 'ISDTYP not in range in reinitialize_current_generator'
        END SELECT

        cg1(g) = lg1(g)
        cg2(g) = lg2(g)
        RETURN

      END SUBROUTINE reinitialize_current_generator

!*********************************************************************

      SUBROUTINE set_all_seeds(iseed1,iseed2)
!----------------------------------------------------------------------
!      SUBROUTINE SETALL(ISEED1,ISEED2)
!               SET ALL random number generators
!     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
!     initial seeds of the other generators are set accordingly, and
!     all generators states are set to these seeds.
!     This is a transcription from Pascal to Fortran of routine
!     Set_Initial_Seed from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     ISEED1 -> First of two integer seeds
!                                   INTEGER ISEED1
!     ISEED2 -> Second of two integer seeds
!                                   INTEGER ISEED1
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: iseed1, iseed2
! ..
! .. Local Scalars ..
        INTEGER :: g, ocgn
! ..
! .. Executable Statements ..

!     TELL THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
!      HAS BEEN CALLED.
        seeds_set = .TRUE.
        CALL get_current_generator(ocgn)
        ig1(1) = iseed1
        ig2(1) = iseed2
        CALL reinitialize_current_generator(-1)
        DO g = 2, n_generators
          ig1(g) = multiply_modulo(a1vw,ig1(g-1),m1)
          ig2(g) = multiply_modulo(a2vw,ig2(g-1),m2)
          CALL set_current_generator(g)
          CALL reinitialize_current_generator(-1)
        END DO
        CALL set_current_generator(ocgn)
        RETURN

      END SUBROUTINE set_all_seeds

!*********************************************************************

      SUBROUTINE set_antithetic(qvalue)
!----------------------------------------------------------------------
!      SUBROUTINE SETANT(QVALUE)
!               SET ANTithetic
!     Sets whether the current generator produces antithetic values.  If
!     X   is  the value  normally returned  from  a uniform [0,1] random
!     number generator then 1  - X is the antithetic  value. If X is the
!     value  normally  returned  from a   uniform  [0,N]  random  number
!     generator then N - 1 - X is the antithetic value.
!     All generators are initialized to NOT generate antithetic values.
!     This is a transcription from Pascal to Fortran of routine
!     Set_Antithetic from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     QVALUE -> .TRUE. if generator G is to generating antithetic
!                    values, otherwise .FALSE.
!                                   LOGICAL QVALUE
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        LOGICAL, INTENT (IN) :: qvalue
! ..
! .. Local Scalars ..
        INTEGER :: g
! ..
! .. Executable Statements ..

        CALL get_current_generator(g)
        qanti(g) = qvalue
        RETURN

      END SUBROUTINE set_antithetic

!*********************************************************************

      SUBROUTINE set_current_generator(g)
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: g
! ..
! .. Executable Statements ..

        IF ((g<1) .OR. (g>n_generators)) THEN
          WRITE (*,90000) g, n_generators
90000     FORMAT ( &
            ' Generator number out of range in set_current_generator'/ &
            ' Generator number is: ',I0/' Allowable range is 1 to : ',I0)
          STOP 'Generator number out of range in set_current_generator'

        ELSE
          current_generator = g
          generator_set = .TRUE.
        END IF

        RETURN

      END SUBROUTINE set_current_generator

!*********************************************************************

      SUBROUTINE set_current_seed(iseed1,iseed2)
!----------------------------------------------------------------------
!               SET S-ee-D of current generator
!     Resets the initial  seed of  the current  generator to  ISEED1 and
!     ISEED2. The seeds of the other generators remain unchanged.
!     This is a transcription from Pascal to Fortran of routine
!     Set_Seed from the paper
!     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
!     with Splitting Facilities." ACM Transactions on Mathematical
!     Software, 17:98-111 (1991)
!                              Arguments
!     ISEED1 -> First integer seed
!                                   INTEGER ISEED1
!     ISEED2 -> Second integer seed
!                                   INTEGER ISEED1
!----------------------------------------------------------------------
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: iseed1, iseed2
! ..
! .. Local Scalars ..
        INTEGER :: g
! ..
! .. Executable Statements ..

        IF ( .NOT. seeds_set) CALL set_default
        CALL get_current_generator(g)
        ig1(g) = iseed1
        ig2(g) = iseed2
        CALL reinitialize_current_generator(-1)
        RETURN

      END SUBROUTINE set_current_seed

!*********************************************************************

      SUBROUTINE set_default
!!! Set base generator to default values and issue warning message
! .. Implicit None Statement ..
        IMPLICIT NONE
! ..
! .. Parameters ..
        CHARACTER (*), PARAMETER :: message_format = '(/5X,''WARNING: &
          &The seeds for the random number generator have been''/5X,''set &
          & to a default  value.  All  random numbers  generated &
          &in''/5X,''different runs of this program that use &
          &this default setting''/5X,''will be the same.''/)'
! ..
! .. Executable Statements ..

        IF ( .NOT. seeds_set) THEN
          CALL set_all_seeds(default_seed1,default_seed2)
          WRITE (*,message_format)
        END IF

        IF ( .NOT. generator_set) CALL set_current_generator(1)

        RETURN

      END SUBROUTINE set_default

!*********************************************************************

    END MODULE ecuyer_cote_mod
