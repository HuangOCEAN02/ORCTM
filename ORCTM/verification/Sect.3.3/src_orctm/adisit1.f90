      SUBROUTINE ADISIT1(TH,SH,PA)
      USE MO_PARAM3
      REAL TH,SH,PA,TPO
!***********************************************************
!     TRANSFORMATION FROM POTENTIAL TO IN SITU DENSITY
!     ACCORDING BRYDEN DSR 20, 401-408 (GILL P.602)
!      WHICH GIVES THE INVERSE TRANSFORMATION
!     FOR AN APPROXIMATE VALUE, ALL TERMS LINEAR IN T ARE TAKEN
!    AFTER THAT ONE NEWTON STEP
!   FOR THE CHECK VALUE 8.4678516     THE ACCURACY IS 0.2 MIKROKELVIN
!CC***CHANGE BY SSD***
! PA SHOULD NOT BE DIVIDED BY 100 IN THE PRESENT VERSION
       PR=PA
!
! USE MODULE MO_PARAM3
!OtB      A1=3.6504 E-4
!OtB      A2=8.3198 E-5
!OtB      A3=5.4065E-7
!OtB      A4=4.0274 E-9
!OtB      B1=1.7439E-5
!OtB      B2=2.9778E-7
!OtB      C1=8.9309E-7
!OtB      C2=3.1628 E-8
!OtB      C3=2.1987E-10
!OtB      D=4.1057 E-9
!OtB      E1=1.6056 E-10
!OtB      E2=5.0484 E-12
!
!     TH(1,1)=8.4678516
!     SH(1,1)= 25.
!     PR=1000.
!
      S = SH
      TPO=TH
      TH=(TH+A1*PR+B1*PR*(S-35.)+C1*PR**2                               &
     & -D*PR**2*(S-35.)-E1*PR**3)                                       &
     & /(1.-A2*PR+B2*PR*(S-35.)+C2*PR**2-E2*PR**3)
!
      T=TH
      FNE=T-PR*(A1+A2*T-A3*T**2+A4*T**3)-PR*(S-35.)*(B1-B2*T)           &
     &-PR**2*(C1-C2*T+C3*T**2) +D*(S-35.)*PR**2-PR**3*(-E1+E2*T)        &
     &    -TPO
       FST=1.-PR*(A2-2.*A3*T+3.*A4*T**2) +PR*(S-35.)*B2                 &
     & +PR**2*(C2-2.*C3*T) -E2*PR**3
      TH=TH-FNE/FST
!     IF(I+J.EQ.2) PRINT*,'NEWTON ', FNE,FST
      RETURN
      END
