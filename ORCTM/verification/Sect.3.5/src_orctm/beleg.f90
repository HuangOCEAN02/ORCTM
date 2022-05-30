      SUBROUTINE BELEG
!****************************************************************
!
!     BBBBB   EEEEEE  L       EEEEEE  GGGGG
!     B    B  E       L       E       G
!     BBBBB   EEEEE   L       EEEEE   GGGGG
!     B    B  E       L       E       G    G
!     BBBBB   EEEEEE  LLLLLL  EEEEEE  GGGGGG
!
!
!*****************************************************************
!
!***********************************************************************
!     BELEG : READS AND COMPUTES START VALUES FOR DIFFERENT FIELDS
!             INCLUDING  A.)       ORIENTATION OF GRID POINTS ON EARTH
!                        B.)       DISTANCES BETWEEN POINTS
!                        C.)       LAYER DEPTH CONFIGURATION
!                        D.)       GRID CONFIGURATION
!
!---------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE

      REAL DMAXO, DMINO, DELZ,  RHO
      INTEGER I, J, K, II, JJ
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     LONGITUDE AND LATITUDE OF GRID POINTS
!
!=========LAENGE UND BREITE
!     ODD FIELDS J=1,3,5,...,JTO-1
!
      DMAXO=-1.E20
      DMINO=1.E20
!
      DO 152 J=1,JE
      DO 152 I=1,IE

      II = I + p_joff
      JJ = J + p_joff

      DTDXPO(I,J)=DT/DLXP(I,J)
      DTDYO(I,J)=DT/DLYP(I,J)
!
      DELZ=(DLXP(I,J)+DLYP(I,J))/2
      IF(DELZ.LE.0.) WRITE(IO_STDOUT,*) 'NULL ',II,JJ
!
      DMAXO=MAX(DMAXO,DELZ)
      DMINO=MIN(DMINO,DELZ)
!
152   CONTINUE
!
      CALL global_min(DMINO)
      CALL global_max(DMAXO)
      WRITE(IO_STDOUT,15290) DMINO/1000.,DMAXO/1000.
15290 FORMAT(' RESOLUTION KM  MIN. : ',F10.3,'  MAX. : ',F10.3,' ODD ')

      DO J=1,JE
      DO I=I_start,IE
      DTDXUO(I,J)=DT/DLXU(I,J)
      DPYO(I,J)=DT/DLYU(I,J)
      ENDDO
      ENDDO

      DO J=J_start,JE
      DO I=1,IE
      DTDXUE(I,J)=DT/DLXV(I,J)
      DPYE(I,J)= DT/DLYV(I,J)
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
!     LAYER DEPTH CONFIGURATION  (IS ACTUALLY DECOFTED IN MPGR PROGRAM
!                              BY : DATA DZW/..../ WHICH GIVES THE
!                              LAYER THICKNESSES STARTING FROM THE TOP)
!
!                               TIESTU(K) : DEPTH OF U-POINT IN LAYER K
!                               TIESTW(K) : DEPTH OF W-POINT OF THE
!                                           UPPER BOUNDARY OF LAYER K
      TIESTU(1) = 0.5 * DZW(1)
      TIESTW(1) = 0.0
!
      DO 44 K=1,KE
      DWI(K)       = 1. / DZW(K)
      TIESTW(K+1)  = TIESTW(K) + DZW(K)
44    CONTINUE
!
!-------------------------
!     CALCULATION OF VECTOR POINT DISTANCES DZ(1,..,KEP)
!
      DO 4544 K=2,KE
      TIESTU(K) = 0.5 * ( TIESTW(K+1) + TIESTW(K) )
4544  CONTINUE
      TIESTU(KEP)=9000.
!
      DZ(1) = TIESTU(1)
      DO 4545 K=2,KEP
4545  DZ(K) = TIESTU(K) - TIESTU(K-1)
!
      DO 4504 K=1,KE
      WRITE(IO_STDOUT,4546)K,TIESTW(K),DZW(K)
      WRITE(IO_STDOUT,4547)K,TIESTU(K),DZ (K)
4504  CONTINUE
!
      WRITE(IO_STDOUT,4546)KEP,TIESTW(K),ZERO
      WRITE(IO_STDOUT,4547)KEP,TIESTU(K),DZ (K)
!
4546  FORMAT(' LAYER ',I2,' W-POINT DEPTH ',F6.0,30X,' THICKNESS : ',   &
     &F6.1)
4547  FORMAT(' LAYER ',I2,' U-POINT DEPTH ',15X,F6.0,' DISTANCE  : ',   &
     &F6.1)
!
!--------------------------- END OF LAYER CONFIGURATION ----------------
!
      WRITE(IO_STDOUT,*)IE_G,JE_G,KE
      WRITE(IO_STDOUT,*)' REFERENCE STRATIFICATION : '
!
      DO 701 K=1,KE 
        DI(K)    = 1. / DZ(K)
        PREFF(K) = 0.1026 * TIESTU(K)
        ROREF(K) = RHO ( SAF(K),TAF(K),PREFF(K) )
        WRITE(IO_STDOUT,70199)K,ROREF(K),SAF(K),TAF(K),PREFF(K)
        DO 701 J=1,JE
        DO 701 I=1,IE
          DVO(I,J,K)=1.E-4
          THO(I,J,K)=TAF(K)
          SAO(I,J,K)=SAF(K)
701   CONTINUE

70199 FORMAT(' LAYER ',I3,' RHO : ',F12.4,' S,T,P ',3F10.5)

      DI(KEP)=1./DZ(KE)

      DO 702 J=1,JE
      DO 702 I=1,IE
        DVO(I,J,KEP)=ZERO
        DVO(I,J,1)=ZERO
702   CONTINUE
!
      RETURN
      END
