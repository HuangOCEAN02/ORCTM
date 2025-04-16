      PROGRAM SETUP_CVCTM
      use chpar         

      REAL(dp),DIMENSION(:,:),ALLOCATABLE::depto,gila_ori,giph_ori,depto_smooth
      REAL(dp) :: deptoma,dummy
      INTEGER, PARAMETER :: io_in_topo=80, in_out_topo=81
      INTEGER,ALLOCATABLE,DIMENSION(:,:)::maskin
      
      ALLOCATE(depto(ie,je),gila_ori(ito+2,jto+2),giph_ori(ito+2,jto+2))
     ALLOCATE(maskin(ie,je),depto_smooth(ie,je))
     
      PRINT *, 'IE = ', IE
      PRINT *, 'JE = ', JE
  
      OPEN (io_in_topo,file='CVTM_topo_new.bin',form='BINARY', &
      &      access='SEQUENTIAL',ACTION='READ')
      OPEN (in_out_topo,file='anta_region',form='UNFORMATTED',&
      &      access='SEQUENTIAL')
  
      PRINT *, 'read topography .. '
      PRINT*, ' vor lon '
      READ(io_in_topo)  ((gila_ori(i,j),i=1,ito+2),j=1,jto+2)
  
      PRINT*, ' vor lat'
      READ(io_in_topo) ((giph_ori(i,j),i=1,ito+2),j=1,jto+2)

      PRINT*, ' vor depto '
      READ(io_in_topo) ((depto(i,j),i=1,ie),j=1,je)
      PRINT *, 'read topography .. done'

      gila_ori = gila_ori / REAL(aradtogra, dp)
      giph_ori = giph_ori / REAL(aradtogra, dp)


      do j=1,je
      do i=1,ie
        if (depto(i,j) .EQ. 0._dp) then
          maskin(i,j)=0
        else
          maskin(i,j)=1
        endif
      enddo
      enddo
      
     ! The Geba Estuary area
      DO I=1,100
        call smooth(depto,maskin,depto_smooth,537,ie-1,396,je-1)
        depto = depto_smooth
      ENDDO

! modify shallow and deep and coast topo start 

      DO j = 1,je
        DO i = 1,ie
           depto(i,j) = MAX(0.,depto(i,j))
           
        IF (depto(i,j) > 0. .AND. depto(i,j) < lev1st) THEN
            depto(i,j) = lev1st
        ENDIF

        IF (depto(i,j) > max_dep) THEN
            depto(i,j) = max_dep
        ENDIF
             
       ENDDO
      ENDDO
    
 ! modify west boundary topo
      DO J=1,JE
         if (depto(1,j)>0.0) then
          if (depto(2,j)>0.0) then
              depto(1,j)=depto(2,j)
             if(depto(3,j)<1.0) depto(3,j)=depto(2,j)
          else
             depto(1,j)=0.0
          endif
         else
          depto(2,j)=0.0
         endif
      ENDDO
 ! modify east boundary topo
      DO J=1,JE
         if (depto(IE,j)>0.0) then
          if (depto(IE-1,j)>0.0) then
              depto(IE,j)=depto(IE-1,j)
             if(depto(IE-2,j)<1.0) depto(IE-2,j)=depto(IE-1,j)
          else
             depto(IE,j)=0.0
          endif
         else
          depto(IE-1,j)=0.0
         endif
      ENDDO
        
 ! modify north boundary topo
      DO I=1,IE
         if (depto(i,1)>1.0) then
          if (depto(i,2)>1.0) then
              depto(i,1)=depto(i,2)
             if(depto(i,3)<1.0) depto(i,3)=depto(i,2)
          else
             depto(i,1)=0.0
          endif
         else
          depto(i,2)=0.0
         endif
      ENDDO
  
 ! modify south boundary topo
      DO I=1,IE
         if (depto(i,JE)>1.0) then
          if (depto(i,JE-1)>1.0) then
              depto(i,JE)=depto(i,JE-1)
              if(depto(i,JE-2)<1.0) depto(i,JE-2)=depto(i,JE-1)
          else
              depto(i,JE)=0.0
          endif
         else
          depto(i,JE-1)=0.0
         endif
      ENDDO

          
      WRITE (in_out_topo) 0_i8,  54_i8, -99_i8, INT((ito+2)*(jto+2),i8)
      WRITE (in_out_topo) ((gila_ori(i,j),i=1,ito+2),j=1,jto+2)

      WRITE (in_out_topo) 0_i8,  55_i8, -99_i8, INT((ito+2)*(jto+2),i8)
      WRITE (in_out_topo) ((giph_ori(i,j),i=1,ito+2),j=1,jto+2)

      WRITE (in_out_topo) 0_i8, 507_i8, -99_i8, INT(ie*je,i8)
      WRITE (in_out_topo) ((depto(i,j),i=1,ie),j=1,je)

      deptoma=0.0
      DO j=1,je
         DO i=1,ie
          deptoma = MAX(depto(i,j),deptoma)
         ENDDO
      ENDDO      
      WRITE (6,*) 'DEPTOMA: ', deptoma

      END PROGRAM SETUP_CVCTM

     subroutine smooth(depto_reg,maskin,depto_pro,i1,i2,j1,j2)
      ! lon(lat) should be set in integal
      ! do not support lon exact location
      
      use chpar
      real(dp)::depto_reg(ie,je),depto_pro(ie,je)
      integer::i,j,i1,i2,j1,j2,ia,ib,ic,ja,jb,jc,n,maskin(ie,je)
          
      depto_pro=depto_reg
     
      PRINT*, 'Smooth NOW!!!'
 
      do j=j1,j2
          do i=i1,i2
            if ((j.eq.1).or.(j.eq.je)) goto 30
              ib=i
              jb=j
              call cal_abc(ia,ib,ic,ja,jb,jc)
              if (maskin(ib,jb) .eq. 1) then
                  n=maskin(ia,jb)+maskin(ic,jb)+maskin(ib,ja)+maskin(ib,jc) &
                   +maskin(ia,ja)+maskin(ic,ja)+maskin(ia,jc)+maskin(ic,jc) &
                   +maskin(ib,jb)
                  if (n < 9) then
                      depto_pro(ib,jb)=(depto_reg(ib,jb)*maskin(ib,jb)+depto_reg(ia,jb)*maskin(ia,jb)&
                                     +depto_reg(ic,jb)*maskin(ic,jb)+depto_reg(ib,ja)*maskin(ib,ja) &
                                     +depto_reg(ib,jc)*maskin(ib,jc)+depto_reg(ia,ja)*maskin(ia,ja) &
                                     +depto_reg(ic,jc)*maskin(ic,jc)+depto_reg(ic,ja)*maskin(ic,ja) &
                                     +depto_reg(ia,jc)*maskin(ia,jc))/n
                  else
                      depto_pro(ib,jb)=depto_reg(ib,jb)/2. &
                   +(depto_reg(ia,jb)+depto_reg(ic,jb)+depto_reg(ib,ja)+depto_reg(ib,jc))/12.&
                   +(depto_reg(ia,ja)+depto_reg(ic,jc)+depto_reg(ic,ja)+depto_reg(ia,jc))/24.
                  endif
              endif
      30    continue
          enddo
      enddo
      PRINT*, 'Smooth END !!!'      
      return
      end subroutine smooth
      
      
       
      subroutine cal_abc(ia,ib,ic,ja,jb,jc)

      integer::ia,ib,ic,ja,jb,jc
      ja=jb-1
      jc=jb+1
      
      ia=ib-1
      ic=ib+1
      
      return
      end subroutine cal_abc
