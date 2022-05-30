      program DoubleRidge_case
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
      INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
      INTEGER, PARAMETER :: spongew = 30
      INTEGER, PARAMETER :: spongee = 30      

      real(dp),dimension(:,:),allocatable::                             &
     & depto,gila,giph,DLXP,DLXU,DLXV,DLYP,DLYU,DLYV,FTWOU,FTWOV,       &
     & UHOPE,VHOPE,TEMOUT,PRESSOUT,TDEWOUT,SWRADOUT,DWLWOUT,PRATEOUT,   &
     & U10OUT,RIV
      real(dp),dimension(:,:,:),allocatable::TEM,SAL
      real(dp),allocatable::TEM_1D(:),SAL_1D(:)
      real(dp),allocatable::xx(:)
      real(dp)::J
      real(dp),parameter::ddx=200
      integer(i8),parameter::ie=4000,je=5,ito=ie*2,jto=je*2,ke=300
      real(dp)::ZZOUT0(ke),ZZOUT(ke)
      REAL(dp), PARAMETER :: api   = 3.1415926
      REAL(dp), PARAMETER :: trot = 86164._dp
      REAL(dp), PARAMETER :: OMEGA = 2._dp * api/trot
      REAL(dp),dimension(:),allocatable::West_ddx 
      REAL(dp),dimension(:),allocatable::East_ddx
      
      
      dep=3000
      allocate(depto(ie,je),gila(ito+1,jto+1),giph(ito+1,jto+1),        &
     &         DLXP(ie,je),DLXU(0:ie,je),DLXV(ie,0:je),FTWOV(ie,0:je),  &
     &         DLYP(ie,je),DLYU(0:ie,je),DLYV(ie,0:je),FTWOU(0:ie,je),  &
     &    UHOPE(0:ie,je),VHOPE(ie,0:je),TEMOUT(ie,je),PRESSOUT(ie,je),  &
     &  TDEWOUT(ie,je),SWRADOUT(ie,je),DWLWOUT(ie,je),PRATEOUT(ie,je),  &
     &  U10OUT(ie,je),RIV(ie,je) )
      allocate(TEM(ie,je,ke),SAL(ie,je,ke),xx(ie),TEM_1D(ke),SAL_1D(ke))
      allocate(West_ddx(spongew))
      allocate(East_ddx(spongee))

      DATA ZZOUT0 /300*10./ 
! grid ----- anta
      gila=0.0
      giph=0.0
! xx ----- stretch for West Sponge bounday Grid        
      OPEN (88,file='ddxRef_west.bin',form='BINARY', &
        &     access='SEQUENTIAL',ACTION='READ')
      OPEN (99,file='ddxRef_east.bin',form='BINARY', &
        &     access='SEQUENTIAL',ACTION='READ')        
        
      DO K=1,spongew
           READ(88) West_ddx(K) 
           print*,'WEST', K , West_ddx(K) 
      ENDDO      

      DO K=1,spongee
           READ(99) East_ddx(K) 
           print*,'EAST', K , East_ddx(K)            
      ENDDO  
      CLOSE(88)
      CLOSE(99)          
! xx ----- x axis
  
      xx(1) = 0._dp
      DO i = 1,ie-1
         xx(i+1) = xx(i) + ddx
      ENDDO
      
      depto(:,:)=0._dp
      DO I=1,IE
        depto(I,3) = dep - 1300_dp*exp(-((xx(i)-350E+3_dp)/20E+3_dp)**2) &
                     - 2500_dp*exp(-((xx(i)-450E+3_dp)/20E+3_dp)**2)
      ENDDO
      
     !  depto(:,3) = dep
   
     ! write(*,"(10f8.5)"),depto(:,3)
      
      OPEN (12,file='anta',form='UNFORMATTED',access='SEQUENTIAL')
      WRITE (12) int(0,i8), int(54,i8), int(-99,i8), (ito+1)*(jto+1)
      WRITE (12) gila
      WRITE (12) int(0,i8), int(55,i8), int(-99,i8), (ito+1)*(jto+1)
      WRITE (12) giph
      WRITE (12) int(0,i8), int(507,i8), int(-99,i8), ie*je
      WRITE (12) depto
      CLOSE (12)

! grid ----- arcgri
      DO i = 1,ie
      	if(i<=spongew) then
         DLXP(i,:)=West_ddx(i)
        ! DLYP(i,:)=West_ddx(i)
        elseif (i>=ie-spongee+1)then
             DLXP(i,:)=East_ddx(i-ie+spongee)
           !  DLYP(i,:)=East_ddx(i-ie+spongee)
            else
             DLXP(i,:)=ddx
           !  DLYP(i,:)=ddx
        endif
      enddo
      
      DO i = 0,ie
      	if(i<=spongew-1) then
          DLXU(i,:)=West_ddx(i+1)
         ! DLYU(i,:)=West_ddx(i+1)
       elseif (i>=ie-spongee+1) then
            DLXU(i,:)=East_ddx(i-ie+spongee)
          !   DLYU(i,:)=East_ddx(i-ie+spongee)
            else
             DLXU(i,:)=ddx
         !    DLYU(i,:)=ddx
         endif
      ENDDO
      
      DO i = 1,ie
      	if( i <= spongew ) then
         DLXV(i,:)=West_ddx(i)
        ! DLYV(i,:)=West_ddx(i)
          elseif ( i >= ie-spongee+1) then
              DLXV(i,:)=East_ddx(i-ie+spongee)
         !     DLYV(i,:)=East_ddx(i-ie+spongee)
            else
              DLXV(i,:)=ddx
         !     DLYV(i,:)=ddx
            endif
      enddo      
      
!      DLXP(:,:)=ddx
      DLYP(:,:)=ddx
      
!      DLXU(:,:)=ddx
      DLYU(:,:)=ddx
      
!      DLXV(:,:)=ddx
      DLYV(:,:)=ddx
      
      FTWOU=2._dp * omega * sind(20.5)
      FTWOV=2._dp * omega * sind(20.5)

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
      TEMOUT=29.2229
      TDEWOUT=29.2229+273.15
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

      sumtiefe=0.
      do K=1,KE
          zzout(K)=sumtiefe+zzout0(K)/2.
          sumtiefe=sumtiefe+zzout0(K)
      enddo

        OPEN (100,file='TRef.bin',form='BINARY', &
        &     access='SEQUENTIAL',ACTION='READ')
        OPEN (101,file='SRef.bin',form='BINARY', &
        &     access='SEQUENTIAL',ACTION='READ')        
 
        DO K=1,KE
           READ(100) TEM_1D(K) 
           READ(101) SAL_1D(K)         
           TEM(:,:,K) = TEM_1D(K)
           SAL(:,:,K) = SAL_1D(k)
        ENDDO
                         
                            
      OPEN(102,FILE='INITEM',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(105,FILE='INISAL',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      DO K=1,KE
        WRITE(102) int(IF1,i8),int(IF22,i8),INT(ZZOUT(K),i8),(IE*JE)
        WRITE(105) int(IF1,i8),int(IF25,i8),INT(ZZOUT(K),i8),(IE*JE)
        WRITE(102) TEM(:,:,k)
        WRITE(105) SAL(:,:,k)
      ENDDO



      end program DoubleRidge_case
