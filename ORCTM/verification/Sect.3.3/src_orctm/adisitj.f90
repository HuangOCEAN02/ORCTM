      SUBROUTINE ADISITJ(TH,SH,PA,J)
! ------------------------------------------------------------------------------
!
!**** *ADISIT*  - TRANSFORMS POTENTIAL TO IN-SITU DENSITY.
!
!     MODIFIED
!     --------
!     O. BOEHRINGER     *DKRZ*                   95
!        - THIS VERSION USES ONLY 44 % OF THE CPU OF THE ORIGINAL HOPC VERSION
!     UWE MIKOLAJEWICZ 2/99
!     ==>ONE-DIMENSIONAL ARRAY, MERGE LOOPS
!
!     METHOD.
!     --------
!     TRANSFORMATION FROM POTENTIAL TO IN SITU DENSITY
!     ACCORDING BRYDEN DSR 20, 401-408 (GILL P.602)
!     WHICH GIVES THE INVERSE TRANSFORMATION
!     FOR AN APPROXIMATE VALUE, ALL TERMS LINEAR IN T ARE TAKEN
!     AFTER THAT ONE NEWTON STEP.
!     FOR THE CHECK VALUE 8.4678516 THE ACCURACY IS 0.2 MIKROKELVIN.
!
!**   INTERFACE.
!     ----------
!     *CALL* *ADISIT(TH,SH,PA)*       CALLED FROM *OCTHER*.
!
!     *COMMON*    *"PARAM1*            - OCEAN GRID DIMENSIONS.
!
!
!     INPUT:
!     -----
!     *TH*        POTENTIAL TEMPERATURE [DEG C]
!     *SH*        SALINITY  [PSU.]
!     *PA*        PRESSURE  [PA]
!    
!     OUTPUT:
!     ------
!     *TH*        IN-SITU  TEMPERATURE [DEG C]
!
! ------------------------------------------------------------------------------
!
      USE MO_PARAM1
      USE MO_PARAM3
!
      DIMENSION TH(IE,JE),SH(IE,JE)
!
      PR=PA
!
!  CHECK VALUES
!     TH(1)=8.4678516
!     SH(1)= 25.
!     PR=1000.
!
      QC = PR*(A1 + PR*(C1 - E1*PR))
      QV = PR*(B1 - D*PR)
      DC = 1. + PR*(-A2 + PR*(C2 -E2*PR))
      DV = B2*PR
      QNQ  = -PR*(-A3 + PR*C3)
      QN3  = -PR*A4
!
      DO  I=1,IE
!
      QVS = QV*(SH(I,J) - 35.) + QC
      DVS = DV*(SH(I,J) - 35.) + DC
      TPO     = TH(I,J)
      TH(I,J) = (TH(I,J) + QVS)/DVS
      T       = TH(I,J)
      FNE     = - QVS + T*(DVS + T*(QNQ + T*QN3)) - TPO
      FST     = DVS + T*(2*QNQ + 3*QN3*T)
      TH(I,J) =TH(I,J)-FNE/FST
      ENDDO
!
      RETURN
      END
