      SUBROUTINE OCVTRO
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
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
     ! if (have_g_is) then
     !  write(IO_STDOUT,*) ' OCVTRO Start West Barotropic velocity'
     !  write(IO_STDOUT,*) 'USO(I_start:3,:) m/s', USO(I_start:3,:)
     !  endif     
    
     ! if ( have_g_ie ) then          
     !  write(IO_STDOUT,*) ' OCVTRO Start East Barotropic velocity'
     !  write(IO_STDOUT,*) 'USO(IE-2:IE,:) m/s', USO(IE-2:IE,:)
     ! endif  
     
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
!     
       CALL bounds_exch(USO)
       CALL bounds_exch(VSE)
!

      K=1
      DO 7715 J=(J_start+1),JE
      DO 7715 I=(I_start+1),IE
        ZO(I,J)=ZO(I,J)+Z1O(I,J)
7715  CONTINUE
!
      CALL bounds_exch(ZO)
       

      RETURN
      END
