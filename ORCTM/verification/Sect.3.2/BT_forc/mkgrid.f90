      program ISW_tank
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
      INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

      real(dp),dimension(:,:),allocatable::                             &
     & depto,gila,giph,DLXP,DLXU,DLXV,DLYP,DLYU,DLYV,FTWOU,FTWOV,       &
     & UHOPE,VHOPE,TEMOUT,PRESSOUT,TDEWOUT,SWRADOUT,DWLWOUT,PRATEOUT,   &
     & U10OUT,RIV
      real(dp),dimension(:,:),allocatable::TEM,SAL
      real(dp),allocatable::ZZOUT0(:),ZZOUT(:)

      real(dp),parameter::ddx=2.E-3,ddz=1.E-3
      integer(i8),parameter::ie=1004,je=5,ito=ie*2,jto=je*2,ke=100

      salx=20
      dep=ddz*ke

      allocate(depto(ie,je),gila(ito+1,jto+1),giph(ito+1,jto+1),        &
     &         DLXP(ie,je),DLXU(0:ie,je),DLXV(ie,0:je),FTWOV(ie,0:je),  &
     &         DLYP(ie,je),DLYU(0:ie,je),DLYV(ie,0:je),FTWOU(0:ie,je),  &
     &    UHOPE(0:ie,je),VHOPE(ie,0:je),TEMOUT(ie,je),PRESSOUT(ie,je),  &
     &  TDEWOUT(ie,je),SWRADOUT(ie,je),DWLWOUT(ie,je),PRATEOUT(ie,je),  &
     &  U10OUT(ie,je),RIV(ie,je) )
      allocate(TEM(ie,je),SAL(ie,je),ZZOUT0(ke),ZZOUT(ke))

! grid ----- anta
      gila=0.0
      giph=0.0

      DO J=1,JE
       DO I=1,IE
        depto(I,J)=dep
        IF ((J.LT.3).OR.(J.GT.(JE-2))) depto(I,J)=0.0
        IF ((I.LT.3).OR.(I.GT.(IE-2))) depto(I,J)=0.0
       ENDDO
      ENDDO

      !write(*,"(10f8.5)"),depto(:,3)
      OPEN (12,file='anta',form='UNFORMATTED',access='SEQUENTIAL')
      WRITE (12) int(0,i8), int(54,i8), int(-99,i8), (ito+1)*(jto+1)
      WRITE (12) gila
      WRITE (12) int(0,i8), int(55,i8), int(-99,i8), (ito+1)*(jto+1)
      WRITE (12) giph
      WRITE (12) int(0,i8), int(507,i8), int(-99,i8), ie*je
      WRITE (12) depto
      CLOSE (12)

! grid ----- arcgri
      DLXP=ddx
      DLXU=ddx
      DLXV=ddx
      DLYP=ddx
      DLYU=ddx
      DLYV=ddx
      FTWOU=0.0
      FTWOV=0.0

      OPEN(82,FILE='arcgri',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      write(82)int(0,i8),int(507,i8),int(-100,i8),ie*je
      write(82)depto
      write(82)int(0,i8),int(85,i8),int(-100,i8),ie*je
      write(82)dlxp
      write(82)int(0,i8),int(185,i8),int(-100,i8),(ie+1)*je
      write(82)dlxu
      write(82)int(0,i8),int(285,i8),int(-100,i8),ie*(je+1)
      write(82)dlxv
      write(82)int(0,i8),int(86,i8),int(-100,i8),ie*je
      write(82)dlyp
      write(82)int(0,i8),int(186,i8),int(-100,i8),(ie+1)*je
      write(82)dlyu
      write(82)int(0,i8),int(286,i8),int(-100,i8),ie*(je+1)
      write(82)dlyv
      write(82)int(0,i8),int(175,i8),int(-100,i8),(ie+1)*je
      write(82)ftwou
      write(82)int(0,i8),int(176,i8),int(-100,i8),ie*(je+1)
      write(82)ftwov

! forcing ----- GI*
      OPEN(52,FILE='GIWIX',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(53,FILE='GIWIY',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(54,FILE='GITEM',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(55,FILE='GITDEW', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(56,FILE='GIPREC', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(57,FILE='GIDWLW', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(58,FILE='GIU10',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(59,FILE='GISWRAD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(61,FILE='GIPRESS',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(12,FILE='GIRIV',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      UHOPE=0.
      VHOPE=0.
      TEMOUT=10.
      TDEWOUT=10+273.15
      PRATEOUT=0.
      U10OUT=0.
      SWRADOUT=0.
      PRESSOUT=0.
      RIV=0.

      IDAY=0
      IF3=0

      IF2=180
      WRITE(52) int(IDAY,i8),int(IF2,i8),int(IF3,i8),((IE+1)*JE)
      WRITE(52) UHOPE
      IF2=181
      WRITE(53) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*(JE+1))
      WRITE(53) VHOPE
      IF2=167
      WRITE(54) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(54) TEMOUT
      IF2=168
      WRITE(55) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(55) TDEWOUT
      IF2=423
      WRITE(56) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(56) PRATEOUT
      IF2=177
      WRITE(57) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(57) DWLWOUT
      IF2=665
      WRITE(58) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(58) U10OUT
      IF2=276
      WRITE(59) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(59) SWRADOUT
      IF2=151
      WRITE(61) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(61) PRESSOUT
      IF2=999
      WRITE(12) int(IDAY,i8),int(IF2,i8),int(IF3,i8),(IE*JE)
      WRITE(12) RIV

! TS ----- INI*
      IF1=0
      IF22=2
      IF25=5

      do k=1,ke
          zzout0(k)=ddz
      enddo
      sumtiefe=0.
      do K=1,KE
          zzout(K)=sumtiefe+zzout0(K)/2.
          sumtiefe=sumtiefe+zzout0(K)
      enddo

      OPEN(102,FILE='INITEM',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(105,FILE='INISAL',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      DO K=1,KE
        DO I=1,IE
          DO J=1,JE
            TEM(I,J)=10.0
            IF (I.GT.(2+salx)) THEN
              IF (K.GT.(12)) THEN
                SAL(I,J)=35.0
              ELSE
                SAL(I,J)=5.0
              ENDIF
            ELSE
              IF (K.GT.(12+50)) THEN
                SAL(I,J)=35.0
              ELSE
                SAL(I,J)=5.0
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        WRITE(102) int(IF1,i8),int(IF22,i8),INT(ZZOUT(K),i8),(IE*JE)
        WRITE(105) int(IF1,i8),int(IF25,i8),INT(ZZOUT(K),i8),(IE*JE)
        WRITE(102) TEM
        WRITE(105) SAL
      ENDDO



      end program ISW_tank
