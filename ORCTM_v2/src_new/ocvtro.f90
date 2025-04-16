      SUBROUTINE OCVTRO(i_step)
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      USE MO_OBCS_NEST
      
      INTEGER	i_step,I,J,K
      REAL DTH,DTHI
!
      DTHI = 0.5 * DTI
      DTH  = 0.5 * DT
!
      DO J=1,JE
      DO I=I_start,IE
      UZO(I,J)=UZO(I,J)*DEUTIO(I,J)
      USO(I,J)=UZO(I,J)
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      VZE(I,J)=VZE(I,J)*DEUTIE(I,J)
      VSE(I,J)=VZE(I,J)
      ENDDO
      ENDDO
!
#ifndef ETACORRECT
       DO J=2,JE1
       DO I=(I_start+1),IE1
       USO(I,J)=UZO(I,J)                                                &
     & +rlsa1*G*STABN*DTDXUO(I,J)*(Z1O(I,J)-Z1O(I+1,J))                 &
     & *AMSUO(I,J,1)
       ENDDO
       ENDDO
  
       DO J=(J_start+1),JE1
       DO I=2,IE1
       VSE(I,J)=VZE(I,J)                                                &
     & +rlsa1*G*STABN*DPYE(I,J)*(Z1O(I,J+1)-Z1O(I,J))                   &
     & *AMSUE(I,J,1)
       ENDDO
       ENDDO

! TST grid should be exclued
#ifdef OBCS_UV_FLATHER 
#ifndef OBCS_TST_NEST
       DO J = 2, JE1

         if ( West .AND. have_g_is) then
           USO(I_start+1,J) = UZO(I_start+1,J)
         endif
         
         if ( East .AND. have_g_ie) then
           USO(IE1,J) = UZO(IE1,J)
         endif
         
       ENDDO

       DO I = 2, IE1

         if ( North .AND. have_g_js) then
           VSE(I,J_start+1) = VZE(I,J_start+1)
         endif
         
         if ( South .AND. have_g_je) then
           VSE(I,JE1) = VZE(I,JE1)
         endif

       ENDDO
#endif
#endif
!     
       CALL bounds_exch(USO)
       CALL bounds_exch(VSE)
!
#endif /*ETACORRECT*/

!  OBCS FOR BT TIDAL LEVEL
!  Updated by H.H on 24.04.2023 in Kiel Germany
!  Open boundary values upgradation
!  Radiation Conditions for Z1O USO VSE 

#ifdef OBCS_ETA_FLATHER 
!  Radiation Conditions for Z1O
       CALL OBCS_CHAPMAN
#endif

     ! IF (have_g_is .AND. WEST) THEN
     !   write(IO_STDOUT,*) 'OCVTRO West'
     !   write(IO_STDOUT,*) 'Z1O(1,2:10)',Z1O(1,2:10)
     !   write(IO_STDOUT,*) 'ZO(1,2:10)',ZO(1,2:10)
     !   write(IO_STDOUT,*) ' The second grid'
     !   WRITE(IO_STDOUT,*)  'Z1O(2,:)', Z1O(2,2:10)
     !   WRITE(IO_STDOUT,*)  'ZO(2,:)', ZO(2,2:10) 
     !   WRITE(IO_STDOUT,*)   'OCVTRO END'
     ! ENDIF

     ! IF (have_g_ie) THEN
     !   write(IO_STDOUT,*) 'TRONEU END East'
     !   write(IO_STDOUT,*) 'Z1O(IE-2:IE,:)',Z1O(IE-2:IE,:)
     ! ENDIF
     
      K=1
      DO 7715 J=(J_start+1),JE
      DO 7715 I=(I_start+1),IE
        ZO(I,J)=ZO(I,J)+Z1O(I,J)
7715  CONTINUE
!
      CALL bounds_exch(ZO)
      
      
!      if (have_g_is) then
!      	WRITE(IO_STDOUT,*) ' OBCS_CHAPMAN West  ', ZO(1,:)
!      endif
      
!      if (have_g_js) then
!      	WRITE(IO_STDOUT,*) ' OBCS_CHAPMAN North  ', ZO(:,1)
!      endif
      
!      if (have_g_je) then
!      	WRITE(IO_STDOUT,*) ' OBCS_CHAPMAN South  ', ZO(:,je)
!      endif
      
#ifdef OBCS_UV_FLATHER 

      ZOOLD = ZO  ! upaded ZO for Flather to BTU & BTV
!  Radiation Conditions for USO VSE
#ifdef OBCS_TST_NEST
        CALL OBCS_TST_UVO(i_step)
#else
      CALL OBCS_FLATHER
#endif
      USOOLD = USO    ! upaded BTU & BTV
      VSEOLD = VSE
#endif

      RETURN
      END
