      SUBROUTINE RHO1J(T,S,P,RH,J)
      USE MO_PARAM1
!********************************************************
! ZUSTANDSGLEICHUNG
! UNTERPROGRAMM NACH ADRIAN GILL (ANHANG)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      DIMENSION S(IE,JE),T(IE,JE),RH(IE,JE)
      DATA B0,B1,B2,B3,B4/8.24493E-1,-4.0899E-3,7.6438E-5,              &
     & -8.2467E-7,5.3875E-9/
      DATA C0,C1,C2/-5.72466E-3,1.0227E-4,-1.6546E-6/
      DATA D0/4.8314E-4/
      DATA A0,A1,A2,A3,A4,A5/999.8426,6.793952E-2,                      &
     & -9.095290E-3,1.001685E-4,-1.120083E-6,6.536332E-9/
      DATA F0,F1,F2,F3/54.6746,-0.603459,1.09987E-2,-6.1670E-5/
      DATA G0,G1,G2/7.944E-2,1.6483E-2,-5.3009E-4/
      DATA AI0,AI1,AI2/2.2838E-3,-1.0981E-5,-1.6078E-6/
      DATA AJ0/1.91075E-4/
      DATA AM0,AM1,AM2/-9.9348E-7,2.0816E-8,9.1697E-10/
      DATA E0,E1,E2,E3,E4/19652.21,148.4206,-2.327105,                  &
     & 1.360477E-2,-5.155288E-5/
      DATA H0,H1,H2,H3/3.239908,1.43713E-3,                             &
     & 1.16092E-4,-5.77905E-7/
      DATA AK0,AK1,AK2/8.50935E-5,-6.12293E-6,5.2787E-8/
!
      DO I=1,IE
      S(I,J)=MAX(S(I,J),28.)
      S3H=SQRT(S(I,J)**3)
!
      RH(I,J)=(A0+T(I,J)*(A1+T(I,J)                                     &
     &       *(A2+T(I,J)*(A3+T(I,J)*(A4+T(I,J)*A5))))                   &
     &       +S(I,J)*(B0+T(I,J)*(B1+T(I,J)                              &
     &      *(B2+T(I,J)*(B3+T(I,J)*B4))))+D0*S(I,J)**2                  &
     &      +S3H*(C0+T(I,J)*(C1+C2*T(I,J))) )                           &
     &       /(1.-P/(P*(                                                &
     &   H0+T(I,J)*(H1+T(I,J)*(H2+T(I,J)*H3))                           &
     &  +S(I,J)*(AI0+T(I,J)*(AI1+AI2*T(I,J)))+AJ0*S3H                   &
     &  +(AK0+T(I,J)*(AK1+T(I,J)*AK2)                                   &
     &  +S(I,J)*(AM0+T(I,J)*(AM1+T(I,J)*AM2)))*P)+                      &
     &    E0+T(I,J)*(E1+T(I,J)*(E2+T(I,J)*(E3+T(I,J)*E4)))              &
     &      +S(I,J)*(F0+T(I,J)*(F1+T(I,J)*(F2+T(I,J)*F3)))              &
     &      +S3H*(G0+T(I,J)*(G1+G2*T(I,J)))))
      ENDDO
!
      RETURN
      END
