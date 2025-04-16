      program setup_initialfields
      use chpar

      INTEGER,parameter::barotropic=0
      real(dp),dimension(:,:),allocatable::depto,gila,giph
      real(dp),dimension(:,:,:),allocatable::TEM,SAL
      real(dp),dimension(:)::ZZOUT0(ke)

      real(dp),allocatable::TEM_1D(:),SAL_1D(:)
      real(dp),allocatable::ZZOUT(:)

      ! DATA ZZOUT0 /50*10.,25*20.,20*50.,9*150.,14*200./
      DATA ZZOUT0 /5*5.,6*15.,9*20.,4*50.,5*100.,4*250.,5*500.,2*1000./
      allocate(depto(ie,je),gila(ito,jto),giph(ito,jto),ZZOUT(ke))
      allocate(TEM(ie,je,ke),SAL(ie,je,ke),TEM_1D(ke),SAL_1D(ke))


! TS ----- INI*
      IF1=0
      IF22=2
      IF25=5

      sumtiefe=0.
      do K=1,KE
          zzout(K)=sumtiefe+zzout0(K)/2.
          sumtiefe=sumtiefe+zzout0(K)
          print*, zzout(k)
      enddo


      OPEN(102,FILE='INITEM',ACCESS='SEQUENTIAL',FORM='UNFORMATTED',action='write')
      OPEN(105,FILE='INISAL',ACCESS='SEQUENTIAL',FORM='UNFORMATTED',action='write')

      if (barotropic == 1 ) then
        
        TEM(:,:,:)=10.0_dp
        SAL(:,:,:)=34.0_dp
        

      else

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
                     
      ENDIF

      DO K=1,KE         
         WRITE(102) int(IF1,i8),int(IF22,i8),INT(ZZOUT(K),i8),int(IE*JE,i8)
         WRITE(105) int(IF1,i8),int(IF25,i8),INT(ZZOUT(K),i8),int(IE*JE,i8)
         WRITE(102) TEM(:,:,K)
         WRITE(105) SAL(:,:,K)
      ENDDO
        


      end program setup_initialfields
