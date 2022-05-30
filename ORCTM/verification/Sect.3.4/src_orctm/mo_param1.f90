      MODULE MO_PARAM1

      IMPLICIT NONE

! Global model dimensions

#ifdef VERSIONT43
      INTEGER, PARAMETER :: IE_G=130,JE_G=211
      INTEGER, PARAMETER :: KBBT=340,ILLT=15900
#endif
#ifdef VERSIONGR30
      INTEGER, PARAMETER :: IE_G=700,JE_G=470
      INTEGER, PARAMETER :: KBBT=218,ILLT=8430
#endif
#ifdef VERSIONGR15
      INTEGER, PARAMETER :: IE_G=256,JE_G=220
      INTEGER, PARAMETER :: KBBT=400,ILLT=36504
#endif
#ifdef VERSIONGR09
      INTEGER, PARAMETER :: IE_G=400,JE_G=338
      INTEGER, PARAMETER :: KBBT=700,ILLT=95000
#endif
#ifdef VERSIONGR03
      INTEGER, PARAMETER :: IE_G=1004,JE_G=5
      INTEGER, PARAMETER :: KBBT=63,ILLT=696
#endif
#ifdef VERSIONGIN
      INTEGER, PARAMETER :: IE_G=182,JE_G=84
      INTEGER, PARAMETER :: KBBT=118,ILLT=7430
#endif
#ifdef LEVELS40
      INTEGER, PARAMETER :: KE=150
#else
#ifdef LEVELS30
      INTEGER, PARAMETER :: KE=30
#else
#ifdef LEVELS23
      INTEGER, PARAMETER :: KE=23
#else
      INTEGER, PARAMETER :: KE=20
#endif /*LEVELS23*/
#endif /*LEVELS40*/
#endif /*LEVELS30*/

! Local dimensions (set by domain decompostion)

      INTEGER IE, JE, I_start, J_start

! Derived Parameters (set in set_param1 below)

      INTEGER ITO, JTO, IE1, IE2, IT3, KEP, KE1, JE1, IT4, JE2, JT2,    &
     &        IT1, JEH, IEJTO, IEJEKE, IEJE

! SOLVER

#ifdef SOR
      INTEGER, PARAMETER :: IELIMI=0
#else
      INTEGER, PARAMETER :: IELIMI=1
#endif

! Update the Periodic Boundary Conditions -2021.03.12 by H.H in Qingdao
! -- Regional model set ICYCLI_X=0, global model set ICYCLI_X=1
! -- Cyclic W-E Boundary Condition(East&West are linked together)
! -- Attention!!! The number of CPU in EW-direction must be larger than two
#ifdef CYCLICWE
      INTEGER, PARAMETER :: ICYCLI_X=1
#else
      INTEGER, PARAMETER :: ICYCLI_X=0
#endif
! --Cyclic N-S Boundary Condition(North&South are linked together)
! -- Attention!!! The number of CPU in NS-direction must be larger than two
#ifdef CYCLICNS  
      INTEGER, PARAMETER :: ICYCLI_Y=1
#else
      INTEGER, PARAMETER :: ICYCLI_Y=0
#endif  
! Update the restart parameter Z3DOC -2020.05.22 by H.H in Luoyang
      INTEGER KBB, ILL, IMM, ILT, Z3DOC

      INTEGER, PARAMETER :: NBOX=9

      INTEGER*8 IBLA,ii1,ii2,ii3,ii4,ii5

      CONTAINS

      SUBROUTINE set_param1

! Sets all derived parameters dpeneding on local settings for IE and JE
! Since the times of FORTRAN 66 are long ago, these parameters
! shouldn't be used any more

      ITO=2*IE
      JTO=2*JE
      IE1=IE-1
      IE2=IE-2
      IT3=ITO+3
      KEP=KE+1
      KE1=KE-1
      JE1=JE-1
      IT4=ITO-4
      JE2=JE-2
      JT2=JTO-2
      IT1=ITO-1
      JEH=JE/2
      IEJTO=IE*JTO
      IEJEKE=IE*JE*KE
      IEJE=IE*JE

#ifdef SOR
      KBB=4
      ILL=9
#else
      KBB=KBBT
      ILL=ILLT
#endif
      IMM=2*KBB+1
      ILT=ILL+KBB

      END SUBROUTINE set_param1
      END MODULE MO_PARAM1
