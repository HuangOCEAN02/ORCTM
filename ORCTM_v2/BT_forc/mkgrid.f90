! This is the regional version of Cape Verde Sea Area 
! Make the anta & arcgri
! Written by H.H 2023.04.26
! In order to satify the difference equation,
! the number of gila_ori & giph_ori are ito+2,jto+2
! The final number of the Grid is ie*je
      PROGRAM GRID_CVCTM
      use chpar
      IMPLICIT NONE   
         
      INTEGER :: i, j, mende
      INTEGER(i8) :: ext_hdr(4)
      INTEGER, PARAMETER :: io_in_topo=80, io_in_topo1=83            
      INTEGER, PARAMETER :: io_out_arcgri=82, io_out_anta=81      
      REAL(wp) :: colaee, colane, colase, colass, colasw, colaum, colaup, colavm, &
           colavp, colazz, cophee, cophne, cophse, cophss, cophsw, cophum, &
           cophup, cophvm, cophvp, cophzz, rmindlxp, rmindlyp, silaee, silane, &
           silase, silass, silasw, silaum, silaup, silavm, silavp, silazz, &
           siphee, siphne, siphse, siphsw, siphss, siphum, siphup, siphvm, &
           siphvp, siphzz, xee, xne, xse, xss, xsw, xum, xup, xvm, xvp, &
           xzz, yee, yne, yse, yss, ysw, yum, yup, yvm, yvp, yzz, zee, &
           zne, zse, zss, zsw, zum, zup, zvm, zvp, zzz         
      real(dp),dimension(:,:),allocatable::gila_ori,giph_ori   
      real(dp),dimension(:,:),allocatable::gila,giph,                  &
     & depto,DLXP,DLXU,DLXV,DLYP,DLYU,DLYV,FTWOU,FTWOV


      allocate(gila_ori(ito+2,jto+2),giph_ori(ito+2,jto+2),depto(ie,je))    
      allocate(DLXP(ie,je),DLXU(ie+1,je),DLXV(ie,je+1),FTWOV(ie,je+1),&
     &         DLYP(ie,je),DLYU(ie+1,je),DLYV(ie,je+1),FTWOU(ie+1,je))

      OPEN(io_in_topo, FILE='anta_region', &
      & ACCESS='SEQUENTIAL',FORM='UNFORMATTED',action='read') 
      OPEN(io_out_anta, FILE='anta', &
      & ACCESS='SEQUENTIAL', FORM='UNFORMATTED',action='write')
      OPEN(io_out_arcgri, FILE='arcgri', &
      & ACCESS='SEQUENTIAL', FORM='UNFORMATTED', action='write')
  
      PRINT*, ' vor gila '
      READ(io_in_topo) ext_hdr
      READ(io_in_topo) gila_ori

      PRINT*, ' vor giph'
      READ(io_in_topo) ext_hdr
      READ(io_in_topo) giph_ori

      PRINT*, ' vor depto '
      READ(io_in_topo) ext_hdr
      READ(io_in_topo) depto
      
      mende = 0
      rmindlxp = 10.e9_wp
      rmindlyp = 10.e9_wp
      
      PRINT*, ' ------- dlxp,dlyp'     
! dlxp dlyp 
      DO j = 1, je
        DO i = 1, ie    	 
    	 ! Giph      
           siphup = sin(giph_ori(2 * i + 2, 2 * j + 1))
           cophup = cos(giph_ori(2 * i + 2, 2 * j + 1))
           siphum = sin(giph_ori(2 * i, 2 * j + 1))
           cophum = cos(giph_ori(2 * i, 2 * j + 1))
           
           siphvp = sin(giph_ori(2 * i + 1, 2 * j + 2))
           cophvp = cos(giph_ori(2 * i + 1, 2 * j + 2))
           siphvm = sin(giph_ori(2 * i + 1, 2 * j))
           cophvm = cos(giph_ori(2 * i + 1, 2 * j)) 
       
       ! Gila
           silaup = sin(gila_ori(2 * i + 2, 2 * j + 1))
           colaup = cos(gila_ori(2 * i + 2, 2 * j + 1))
           silaum = sin(gila_ori(2 * i, 2 * j + 1))
           colaum = cos(gila_ori(2 * i, 2 * j + 1))
           
           silavp = sin(gila_ori(2 * i + 1, 2 * j + 2))
           colavp = cos(gila_ori(2 * i + 1, 2 * j + 2))
           silavm = sin(gila_ori(2 * i + 1, 2 * j))
           colavm = cos(gila_ori(2 * i + 1, 2 * j))

           xup = silaup * cophup
           yup = colaup * cophup
           zup = siphup

           xum = silaum * cophum
           yum = colaum * cophum
           zum = siphum

           xvp = silavp * cophvp
           yvp = colavp * cophvp
           zvp = siphvp

           xvm = silavm * cophvm
           yvm = colavm * cophvm
           zvm = siphvm

           dlxp(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xup * xum + yup * yum + zup * zum), 1._wp)))
           
           dlyp(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xvp * xvm + yvp * yvm + zvp * zvm), 1._wp)))
       
           IF (depto(i, j) .GT. 0.5_wp) THEN
      	       mende = mende + 1
               rmindlxp = min(dlxp(i, j), rmindlxp)
               rmindlyp = min(dlyp(i, j), rmindlyp)
           ENDIF
              
        ENDDO
      ENDDO

      PRINT *, 'min dlxp ', rmindlxp
      PRINT *, 'min dlyp ', rmindlyp
      PRINT *, ' sum of depth points ', mende
      

      PRINT*, ' ------- dlxu,dlyu' 
! dlxu dlyu  attention: ee
      DO j = 1, je
        DO i = 1, ie+1 
    	     	 
    	     cophzz = cos(giph_ori(2 * i - 1, 2 * j + 1))
    	     cophee = cos(giph_ori(2 * i + 1, 2 * j + 1))   	 
    	     cophne = cos(giph_ori(2 * i, 2 * j))
    	     cophse = cos(giph_ori(2 * i, 2 * j + 2))
        
           siphzz = sin(giph_ori(2 * i - 1, 2 * j + 1))
    	     siphee = sin(giph_ori(2 * i + 1, 2 * j + 1))
           siphne = sin(giph_ori(2 * i, 2 * j))    	
           siphse = sin(giph_ori(2 * i, 2 * j + 2))
 
           colazz = cos(gila_ori(2 * i - 1, 2 * j + 1))
           colaee = cos(gila_ori(2 * i + 1, 2 * j + 1))
           colane = cos(gila_ori(2 * i, 2 * j))
           colase = cos(gila_ori(2 * i, 2 * j + 2))
       
           silazz = sin(gila_ori(2 * i - 1, 2 * j + 1))
           silaee = SIN(gila_ori(2 * i + 1, 2 * j + 1))
           silane = sin(gila_ori(2 * i, 2 * j))
           silase = sin(gila_ori(2 * i, 2 * j + 2))
       
           xzz = silazz * cophzz
           yzz = colazz * cophzz
           zzz = siphzz

           xee = silaee * cophee
           yee = colaee * cophee
           zee = siphee
 
           xse = silase * cophse
           yse = colase * cophse
           zse = siphse

           xne = silane * cophne
           yne = colane * cophne
           zne = siphne 
       
           dlyu(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xne * xse + yne * yse + zne * zse), 1._wp)))   
         
           dlxu(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xzz * xee + yzz * yee + zzz * zee), 1._wp)))
        
           ftwou(i, j) = 2._wp * omega * sin(giph_ori(2*i , 2*j + 1)) 
        ENDDO
      ENDDO

      PRINT*, 'dlxv,dlyv' 
! dlxv dlyv  attention: ss
      DO j = 1, je+1
        DO i = 1, ie
    	 ! Giph      
           cophzz = cos(giph_ori(2 * i + 1, 2 * j - 1))
           cophsw = cos(giph_ori(2 * i, 2 * j))
           cophss = cos(giph_ori(2 * i + 1, 2 * j + 1))
           cophse = cos(giph_ori(2 * i + 2, 2 * j))

           siphzz = sin(giph_ori(2 * i + 1, 2 * j - 1))
           siphsw = sin(giph_ori(2 * i, 2 * j))
           siphss = sin(giph_ori(2 * i + 1, 2 * j + 1))
           siphse = sin(giph_ori(2 * i + 2, 2 * j))

           colazz = cos(gila_ori(2 * i + 1, 2 * j - 1))
           colasw = cos(gila_ori(2 * i, 2 * j))
           colass = cos(gila_ori(2 * i + 1, 2 * j + 1))
           colase = cos(gila_ori(2 * i + 2, 2 * j))

           silazz = sin(gila_ori(2 * i + 1, 2 * j - 1))
           silasw = sin(gila_ori(2 * i, 2 * j))
           silass = sin(gila_ori(2 * i + 1, 2 * j + 1))
           silase = sin(gila_ori(2 * i + 2, 2 * j))

           xzz = silazz * cophzz
           yzz = colazz * cophzz
           zzz = siphzz

           xss = silass * cophss
           yss = colass * cophss
           zss = siphss

           xsw = silasw * cophsw
           ysw = colasw * cophsw
           zsw = siphsw

           xse = silase * cophse
           yse = colase * cophse
           zse = siphse
     
           dlxv(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xsw * xse + ysw * yse + zsw * zse), 1._wp)))
        
           dlyv(i, j) = MAX(1._wp, radius * &
           &     ACOS(MIN((xzz * xss + yzz * yss + zzz * zss), 1._wp)))
        
           ftwov(i, j) = 2._wp * omega * sin(giph_ori(2*i + 1, 2*j))
      
!      IF (dlxv(i, j) .LT. 1._wp) THEN
!      	  WRITE(6, *)'scheiss dlxv:', i, j, dlxv(i, j), &
!           dlyv(i, j), dlxp(i, j), dlxu(i, j), gila_ori(2*i+1, 2*j+1), giph_ori(2*i+1, 2*j+1)
!      ENDIF
      
        ENDDO
      ENDDO

! Write the arcgri
      ext_hdr(:)=(/0, 193, - 100, ie * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) depto
    
      ext_hdr(:)=(/0, 85, - 100, ie * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlxp

      ext_hdr(:)=(/0, 185, - 100, (ie+1) * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlxu

      ext_hdr(:)=(/0, 188, - 100, ie * (je+1)/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlxv

      ext_hdr(:)=(/0, 86, - 100, ie * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlyp

      ext_hdr(:)=(/0, 186, - 100, (ie+1) * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlyu

      ext_hdr(:)=(/0, 189, - 100, ie * (je+1)/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) dlyv

      ext_hdr(:)=(/0, 132, - 100, (ie+1) * je/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) ftwou

      ext_hdr(:)=(/0, 133, - 100, ie * (je+1)/)
      WRITE(io_out_arcgri) ext_hdr
      WRITE(io_out_arcgri) ftwov

! Upgrade the anta
      WRITE(io_out_anta) 0_i8,  54_i8, -99_i8, INT((ie*2+1)*(je*2+1),i8)
      WRITE(io_out_anta) ((gila_ori(i,j),i=2,ito+1),j=2,jto+1)

      WRITE(io_out_anta) 0_i8,  55_i8, -99_i8, INT((ie*2+1)*(je*2+1),i8)
      WRITE(io_out_anta) ((giph_ori(i,j),i=2,ito+1),j=2,jto+1)

      WRITE(io_out_anta) 0_i8, 507_i8, -99_i8, INT(ie*je,i8)
      WRITE(io_out_anta) ((depto(i,j),i=1,ie),j=1,je)

      end program GRID_CVCTM
