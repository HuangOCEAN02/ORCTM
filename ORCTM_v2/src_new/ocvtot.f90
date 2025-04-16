      SUBROUTINE OCVTOT
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      use mo_obcs_nest
      DIMENSION PSIZ(KEP)

!$OMP PARALLEL PRIVATE(i,j,k)

!$OMP DO
      DO 3130 K=1,KE
!
      DO J=1,JE
      DO I=I_start,IE
      T1O(I,J,K)=UOO(I,J,K)
      UKO(I,J,K)=AMSUO(I,J,K)*(USO(I,J)+UKO(I,J,K))
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      S1O(I,J,K)=VOE(I,J,K)
      VKE(I,J,K)=AMSUE(I,J,K)*(VSE(I,J)+VKE(I,J,K))
      ENDDO
      ENDDO
!
3130  CONTINUE
      CALL bounds_exch(T1O)
      CALL bounds_exch(S1O)
      CALL bounds_exch(UKO)
      CALL bounds_exch(VKE)
!
!$OMP DO
      DO 8001 J=1,JE
      DO 8001 I=1,IE
      WO(I,J,KEP) = ZERO
8001  CONTINUE
!
!$OMP DO
      DO 3170 K=1,KE
!
      DO J=1,JE
      DO I=I_start,IE
      UAKO=UKO(I,J,K)
      UKO(I,J,K)=CONN*UKO(I,J,K)+UOO(I,J,K)*CONO
      UOO(I,J,K)=UAKO
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      VAKE=VKE(I,J,K)
      VKE(I,J,K)=CONN*VKE(I,J,K)+VOE(I,J,K)*CONO
      VOE(I,J,K)=VAKE
      ENDDO
      ENDDO
!
3170  CONTINUE
!
      CALL bounds_exch(UKO)
      CALL bounds_exch(UOO)
      CALL bounds_exch(VKE)
      CALL bounds_exch(VOE)
!
!======================================================================
!
!     B)
!
!     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
!                             HORIZONTAL VELOCITY FIELD
!
!$OMP DO
      DO 7811 K=KEP,1,-1
      DO 7811 J=1,JE
      DO 7811 I=1,IE
      WO(I,J,K) = ZERO
7811  CONTINUE
!
!$OMP DO
      DO 781 K=KE,1,-1
      DO 7815 J=(J_start+1),JE
      DO 7815 I=(I_start+1),IE
!
      WO(I,J,K) = WO(I,J,K+1)                                           &
     & + DTI * WETO(I,J,K) * (                                          &
     & DTDXPO(I,J)   * (   UKO(I-1,J,K) * DDUO(I-1,J,K)*DLYU(I-1,J)     &
     &  - UKO(I,J,K)   * DDUO(I,J,K)*DLYU(I,J)   )/DLYP(I,J)  )         &
     & + DTI*WETO(I,J,K)* (                                             &
     & + DTDYO(I,J) *                                                   &
     & (   VKE(I,J,K)     * DDUE(I,J,K)*DLXV(I,J)                       &
     & - VKE(I,J-1,K)   * DDUE(I,J-1,K)*DLXV(I,J-1)  )/DLXP(I,J)  )
!
7815  CONTINUE
781    CONTINUE
!$OMP END PARALLEL
!
#ifdef OBC_UV_FORC
     ! CALL OBCS_APPLY_WO
#endif 

      CALL bounds_exch(WO)
!
#ifdef NON_HYDROSTATIC
      IF (FIRST_STEP) THEN
        DO 9001 K=1,KEP
        DO 9001 I=1,IE
        DO 9001 J=1,JE  
#ifdef NYDCORRECT
! The 3-D Nonhydrostatic simulation need to be modified
          WNO(I,J,K) = 0.25*WNO(I,J,K) + 0.75*WO(I,J,K)

! This is the bug to be corrected
!        WO(I,J,K) = WNO(I,J,K)  
! ifnot comment, the advection of U&V is very strange !!!
#else
          WO(I,J,K) = WNO(I,J,K)
#endif
 9001 CONTINUE
 
#ifdef OBC_UV_FORC
     CALL OBCS_APPLY_WNO
#endif 
 
      ENDIF
      
      
#endif

      RETURN
      END
