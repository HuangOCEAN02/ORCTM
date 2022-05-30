      SUBROUTINE OCTDIFF_TS
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_UNITS
      USE MO_OCTDIFF
      USE MO_LEVITUS

!
! 
!     COMPUTES DIFFUSION OF TEMPERATURE THO AND SALINITY SAO
!
!     UWE MIKOLAJEWICZ 12/99
!
!     VERTICAL DIFFUSION IMPLICITE
!     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
!                                AH00         AULAPTS
!     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
!                                AHX,AHY      AULX,AULY
!
!     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
!
!     Changes R. Johanni, 2003-11-10:
!     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
!
!     Changes R. Smith, 2004-09-27
!     tracer-independent matrices calculated separately first for use
!     with multiple tracers

!
      implicit none
      INTEGER i, j, k, imal, l
      REAL fakres

      REAL ZZSURF(KE)
      REAL TV_diff, SV_diff, DLXYP,DLXYP1
      REAL TFLUX,TFLUZ,TFLUY,SFLUX,SFLUZ,SFLUY
      REAL THTEN(ie,je,ke),SATEN(ie,je,ke)

      DO K=1,KE
       ZZSURF(K)=0.
      ENDDO    
      ZZSURF(1)=1.

      DO IMAL=1,1

!$OMP PARALLEL PRIVATE(i,j,k,SFLUX,SFLUZ,SFLUY,TFLUX,TFLUZ,TFLUY,FAKRES,TV_diff,SV_diff,DLXYP,DLXYP1)

!put in the tracer values here
!$OMP DO
      DO J=1,JE
       DO I=1,IE
        DO K=1,KE
         T1O(i,j,k)=THO(i,j,k)*vol_term(i,j,k)
         S1O(i,j,k)=SAO(i,j,k)*vol_term(i,j,k)
        ENDDO
       ENDDO
      ENDDO

      IF(AH00.GT.ALMZER)THEN
!
!X DIRECTION - tracer dependent
!$OMP DO
      DO J=2,JE1
       DO I=1,IE1
        DO K=1,KE
         IF(AMSUO(I,J,K).GT.0.)THEN

         TV_diff=(THO(I+1,J,K)-THO(I,J,K))/DLXU(I,J)
         SV_diff=(SAO(I+1,J,K)-SAO(I,J,K))/DLXU(I,J)
         TFLUX=0
         SFLUX=0
#ifdef ISOPYK
         DLXYP=DLXP(I,J)*DLYP(I,J)
         DLXYP1=DLXP(I+1,J)*DLYP(I+1,J)
         TFLUZ=0
         SFLUZ=0
!
!        TRIANGLE LEFT,UPW
!
         TFLUX=TFLUX+XCOEFF_LO(I,J,K)*(TV_diff+SLOPLOx(I,J,K)*          &
     &        (THO(I,J,KO(J,K))-THO(I,J,K))/DZ(K))
         TFLUZ=ZCOEFF_LOx(I,J,K)*(SLOPLOx(I,J,K)*                       &
     &        (THO(I,J,KO(J,K))-THO(I,J,K))*                            &
     &         DLXYP/DZ(K)+1.*DLXYP*TV_diff)
         SFLUX=SFLUX+XCOEFF_LO(I,J,K)*(SV_diff+SLOPLOx(I,J,K)*          &
     &        (SAO(I,J,KO(J,K))-SAO(I,J,K))/DZ(K))
         SFLUZ=ZCOEFF_LOx(I,J,K)*(SLOPLOx(I,J,K)*                       &
     &        (SAO(I,J,KO(J,K))-SAO(I,J,K))*                            &
     &         DLXYP/DZ(K)+1.*DLXYP*SV_diff)

         T1O(I,J,KO(J,K))=T1O(I,J,KO(J,K))-TFLUZ
         S1O(I,J,KO(J,K))=S1O(I,J,KO(J,K))-SFLUZ
         T1O(I,J,K)=T1O(I,J,K)+TFLUZ
         S1O(I,J,K)=S1O(I,J,K)+SFLUZ
!
!        TRIANGLE LEFT DOWN
!
         TFLUX=TFLUX+XCOEFF_LU(I,J,K)*(TV_diff+SLOPLUx(I,J,K)*          &
     &   (THO(I,J,K)-THO(I,J,KU(J,K)))/DZ(KU(J,K)))            
         TFLUZ=ZCOEFF_LUx(I,J,K)*(SLOPLUx(I,J,K)*                       &
     &   (THO(I,J,K)-THO(I,J,KU(J,K)))*                                 &
     &    DLXYP/DZ(KU(J,K))+1.*DLXYP*TV_diff)
         SFLUX=SFLUX+XCOEFF_LU(I,J,K)*(SV_diff+SLOPLUx(I,J,K)*          &
     &   (SAO(I,J,K)-SAO(I,J,KU(J,K)))/DZ(KU(J,K)))            
         SFLUZ=ZCOEFF_LUx(I,J,K)*(SLOPLUx(I,J,K)*                       &
     &   (SAO(I,J,K)-SAO(I,J,KU(J,K)))*                                 &
     &    DLXYP/DZ(KU(J,K))+1.*DLXYP*SV_diff)

         T1O(I,J,K)=T1O(I,J,K)-TFLUZ
         S1O(I,J,K)=S1O(I,J,K)-SFLUZ
         T1O(I,J,KU(J,K))=T1O(I,J,KU(J,K))+TFLUZ
         S1O(I,J,KU(J,K))=S1O(I,J,KU(J,K))+SFLUZ
!
!        TRIANGLE RIGHT,UPW
!
         TFLUX=TFLUX+XCOEFF_RO(I,J,K)*(TV_diff+SLOPROx(I,J,K)*          &
     &   (THO(I+1,J,KO(J,K))-THO(I+1,J,K))/DZ(K))
         TFLUZ=ZCOEFF_ROx(I,J,K)*(SLOPROx(I,J,K)*                       &
     &   (THO(I+1,J,KO(J,K))-THO(I+1,J,K))*                             &
     &    DLXYP1/DZ(K)+1.*DLXYP1*TV_diff)
         SFLUX=SFLUX+XCOEFF_RO(I,J,K)*(SV_diff+SLOPROx(I,J,K)*          &
     &   (SAO(I+1,J,KO(J,K))-SAO(I+1,J,K))/DZ(K))
         SFLUZ=ZCOEFF_ROx(I,J,K)*(SLOPROx(I,J,K)*                       &
     &   (SAO(I+1,J,KO(J,K))-SAO(I+1,J,K))*                             &
     &    DLXYP1/DZ(K)+1.*DLXYP1*SV_diff)

         T1O(I+1,J,KO(j,K))=T1O(I+1,J,KO(j,k))-TFLUZ
         S1O(I+1,J,KO(j,K))=S1O(I+1,J,KO(j,k))-SFLUZ
         T1O(I+1,J,K)=T1O(I+1,J,K)+TFLUZ
         S1O(I+1,J,K)=S1O(I+1,J,K)+SFLUZ
!
!        TRIANGLE RIGHT DOWN
!
         TFLUX=TFLUX+XCOEFF_RU(I,J,K)*(TV_diff+SLOPRUx(I,J,K)*          &
     &   (THO(I+1,J,K)-THO(I+1,J,KU(J,K)))/DZ(KU(J,K)))
         TFLUZ=ZCOEFF_RUx(I,J,K)*(SLOPRUx(I,J,K)*                       &
     &   (THO(I+1,J,K)-THO(I+1,J,KU(J,K)))*                             &
     &    DLXYP1/DZ(KU(J,K))+1.*DLXYP1*TV_diff)
         SFLUX=SFLUX+XCOEFF_RU(I,J,K)*(SV_diff+SLOPRUx(I,J,K)*          &
     &   (SAO(I+1,J,K)-SAO(I+1,J,KU(J,K)))/DZ(KU(J,K)))
         SFLUZ=ZCOEFF_RUx(I,J,K)*(SLOPRUx(I,J,K)*                       &
     &   (SAO(I+1,J,K)-SAO(I+1,J,KU(J,K)))*                             &
     &    DLXYP1/DZ(KU(J,K))+1.*DLXYP1*SV_diff)

         T1O(I+1,J,K)=T1O(I+1,J,K)-TFLUZ
         S1O(I+1,J,K)=S1O(I+1,J,K)-SFLUZ
         T1O(I+1,J,KU(j,K))=T1O(I+1,J,KU(j,k))+TFLUZ
         S1O(I+1,J,KU(j,K))=S1O(I+1,J,KU(j,k))+SFLUZ
#else
         TFLUX=TV_diff*XFLUX(i,j,k)
         SFLUX=SV_diff*XFLUX(i,j,k)
#endif
         T1O(I,J,K)=T1O(I,J,K)+TFLUX
         S1O(I,J,K)=S1O(I,J,K)+SFLUX
         T1O(I+1,J,K)=T1O(I+1,J,K)-TFLUX
         S1O(I+1,J,K)=S1O(I+1,J,K)-SFLUX

         ENDIF
        ENDDO
       ENDDO
      ENDDO

!Y DIRECTION - tracer dependent values
!$OMP DO
      DO J=1,JE1
       DO I=2,IE1
        DO K=1,KE
         IF(AMSUE(I,J,K).GT.0.)THEN

         TV_diff=(THO(I,J+1,K)-THO(I,J,K))/DLYV(I,J)
         SV_diff=(SAO(I,J+1,K)-SAO(I,J,K))/DLYV(I,J)
         TFLUY=0
         SFLUY=0
#ifdef ISOPYK
         DLXYP=DLXP(I,J)*DLYP(I,J)
         DLXYP1=DLXP(I,J+1)*DLYP(I,J+1)
         TFLUZ=0
         SFLUZ=0
!
!        TRIANGLE LEFT,UPW
!
         TFLUY=TFLUY+YCOEFF_LO(I,J,K)*(TV_diff+SLOPLOy(I,J,K)*          &
     &      (THO(I,J,KO(J,K))-THO(I,J,K))/DZ(K))
         TFLUZ=ZCOEFF_LOy(I,J,K)*(SLOPLOy(I,J,K)*                       &
     &      (THO(I,J,KO(J,K))-THO(I,J,K))*                              &
     &       DLXYP/DZ(K)+1.*DLXYP*TV_diff)
         SFLUY=SFLUY+YCOEFF_LO(I,J,K)*(SV_diff+SLOPLOy(I,J,K)*          &
     &      (SAO(I,J,KO(J,K))-SAO(I,J,K))/DZ(K))
         SFLUZ=ZCOEFF_LOy(I,J,K)*(SLOPLOy(I,J,K)*                       &
     &      (SAO(I,J,KO(J,K))-SAO(I,J,K))*                              &
     &       DLXYP/DZ(K)+1.*DLXYP*SV_diff)

         T1O(I,J,KO(J,K))=T1O(I,J,KO(J,K))-TFLUZ
         S1O(I,J,KO(J,K))=S1O(I,J,KO(J,K))-SFLUZ
         T1O(I,J,K)=T1O(I,J,K)+TFLUZ
         S1O(I,J,K)=S1O(I,J,K)+SFLUZ
!
!        TRIANGLE LEFT DOWN
!
         TFLUY=TFLUY+YCOEFF_LU(I,J,K)*(TV_diff+SLOPLUy(I,J,K)*          &
     &    (THO(I,J,K)-THO(I,J,KU(J,K)))/DZ(KU(J,K)))
         TFLUZ=ZCOEFF_LUy(I,J,K)*(SLOPLUy(I,J,K)*                       &
     &    (THO(I,J,K)-THO(I,J,KU(J,K)))*                                &
     &     DLXYP/DZ(KU(J,K))+1.*DLXYP*TV_diff)
         SFLUY=SFLUY+YCOEFF_LU(I,J,K)*(SV_diff+SLOPLUy(I,J,K)*          &
     &    (SAO(I,J,K)-SAO(I,J,KU(J,K)))/DZ(KU(J,K)))
         SFLUZ=ZCOEFF_LUy(I,J,K)*(SLOPLUy(I,J,K)*                       &
     &    (SAO(I,J,K)-SAO(I,J,KU(J,K)))*                                &
     &     DLXYP/DZ(KU(J,K))+1.*DLXYP*SV_diff)

         T1O(I,J,K)=T1O(I,J,K)-TFLUZ
         S1O(I,J,K)=S1O(I,J,K)-SFLUZ
         T1O(I,J,KU(J,K))=T1O(I,J,KU(J,K))+TFLUZ
         S1O(I,J,KU(J,K))=S1O(I,J,KU(J,K))+SFLUZ
!
!        TRIANGLE RIGHT,UPW
!
         TFLUY=TFLUY+YCOEFF_RO(I,J,K)*(TV_diff+SLOPROy(I,J,K)*          &
     &       (THO(I,J+1,KO(J,K))-THO(I,J+1,K))/DZ(K))
         TFLUZ=ZCOEFF_ROy(I,J,K)*(SLOPROy(I,J,K)*                       &
     &       (THO(I,J+1,KO(J,K))-THO(I,J+1,K))*                         &
     &        DLXYP1/DZ(K)+1.*DLXYP1*TV_diff)
         SFLUY=SFLUY+YCOEFF_RO(I,J,K)*(SV_diff+SLOPROy(I,J,K)*          &
     &       (SAO(I,J+1,KO(J,K))-SAO(I,J+1,K))/DZ(K))
         SFLUZ=ZCOEFF_ROy(I,J,K)*(SLOPROy(I,J,K)*                       &
     &       (SAO(I,J+1,KO(J,K))-SAO(I,J+1,K))*                         &
     &        DLXYP1/DZ(K)+1.*DLXYP1*SV_diff)

         T1O(I,J+1,KO(j,K))=T1O(I,J+1,KO(j,k))-TFLUZ
         S1O(I,J+1,KO(j,K))=S1O(I,J+1,KO(j,k))-SFLUZ
         T1O(I,J+1,K)=T1O(I,J+1,K)+TFLUZ
         S1O(I,J+1,K)=S1O(I,J+1,K)+SFLUZ
!
!        TRIANGLE RIGHT DOWN
!
         TFLUY=TFLUY+YCOEFF_RU(I,J,K)*(TV_diff+SLOPRUy(I,J,K)*          &
     &   (THO(I,J+1,K)-THO(I,J+1,KU(J,K)))/DZ(KU(J,K)))
         TFLUZ=ZCOEFF_RUy(I,J,K)*(SLOPRUy(I,J,K)*                       &
     &   (THO(I,J+1,K)-THO(I,J+1,KU(J,K)))*                             &
     &    DLXYP1/DZ(KU(J,K))+1.*DLXYP1*TV_diff)
         SFLUY=SFLUY+YCOEFF_RU(I,J,K)*(SV_diff+SLOPRUy(I,J,K)*          &
     &   (SAO(I,J+1,K)-SAO(I,J+1,KU(J,K)))/DZ(KU(J,K)))
         SFLUZ=ZCOEFF_RUy(I,J,K)*(SLOPRUy(I,J,K)*                       &
     &   (SAO(I,J+1,K)-SAO(I,J+1,KU(J,K)))*                             &
     &    DLXYP1/DZ(KU(J,K))+1.*DLXYP1*SV_diff)

         T1O(I,J+1,K)=T1O(I,J+1,K)-TFLUZ
         S1O(I,J+1,K)=S1O(I,J+1,K)-SFLUZ
         T1O(I,J+1,KU(j,K))=T1O(I,J+1,KU(j,k))+TFLUZ
         S1O(I,J+1,KU(j,K))=S1O(I,J+1,KU(j,k))+SFLUZ
#else      
         TFLUY=YFLUX(i,j,k)*TV_diff
         SFLUY=YFLUX(i,j,k)*SV_diff
#endif
         T1O(I,J,K)=T1O(I,J,K)+TFLUY
         S1O(I,J,K)=S1O(I,J,K)+SFLUY
         T1O(I,J+1,K)=T1O(I,J+1,K)-TFLUY
         S1O(I,J+1,K)=S1O(I,J+1,K)-SFLUY

         ENDIF
        ENDDO
       ENDDO
      ENDDO   

!$OMP SINGLE
      CALL bounds_exch(T1O)
      CALL bounds_exch(S1O)
!$OMP END SINGLE

!$OMP DO
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
          THTEN(i,j,k)=T1O(i,j,k)/vol_term(i,j,k)-THO(i,j,k)
          SATEN(i,j,k)=S1O(i,j,k)/vol_term(i,j,k)-SAO(i,j,k)
        ENDDO
       ENDDO
      ENDDO
!
      ELSE ! AH00.GT.ALMZER
!$OMP DO
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
          THTEN(i,j,k)=0.
          SATEN(i,j,k)=0.
        ENDDO
       ENDDO
      ENDDO
      ENDIF ! AH00.GT.ALMZER

! VERTICAL DIFFUSION -tracer dependent
!$OMP DO
      DO J=2,JE1
      DO I=2,IE1
!
        DO K=1,KE
          T1O(I,J,K)=THO(i,j,k)
          S1O(I,J,K)=SAO(i,j,k)
        ENDDO
!
!     IMPLICIT VERTICAL DIFFUSION
!
        DO K=2,KE
          T1O(I,J,K)=T1O(I,J,K)-TRIDSY(I,J,K-1,1)*T1O(I,J,K-1)
          S1O(I,J,K)=S1O(I,J,K)-TRIDSY(I,J,K-1,1)*S1O(I,J,K-1)
        ENDDO
!
        K=KE
        THO(I,J,K)=T1O(I,J,K)/TRIDSY(I,J,K,2)
        SAO(I,J,K)=S1O(I,J,K)/TRIDSY(I,J,K,2)
!
        DO K=1,KE1
          L=KE-K
          THO(I,J,L)=( T1O(I,J,L)-TRIDSY(I,J,L,3)*THO(I,J,L+1) )        &
     &                / TRIDSY(I,J,L,2)
          SAO(I,J,L)=( S1O(I,J,L)-TRIDSY(I,J,L,3)*SAO(I,J,L+1) )        &
     &                / TRIDSY(I,J,L,2)
        ENDDO

      ENDDO
      ENDDO

!$OMP SINGLE
       CALL bounds_exch(THO)
       CALL bounds_exch(SAO)
!$OMP END SINGLE

!$OMP DO
       DO k=1,ke
        DO j=2,je1
         DO i=2,ie1
          THO(i,j,k)=THO(i,j,k)+THTEN(I,J,K)
          SAO(i,j,k)=SAO(i,j,k)+SATEN(I,J,K)
         ENDDO
        ENDDO
       ENDDO

!
! BIHARMONIC T,S DIFFUSION - tracer dependent
!
       IF(AULAPTS.GT.ALMZER)THEN
!$OMP DO
      do k=1,ke
      do j=2,je1
      do i=2,ie1
      T1O(i,j,k)=weto(i,j,k)*(                                          &
     & (weto(i-1,j,k)*(THO(i-1,j,k)-THO(i,j,k))/dlxu(i-1,j)             &
     &+weto(i+1,j,k)*(THO(i+1,j,k)-THO(i,j,k))/dlxu(i,j))/dlxp(i,j)     &
     &+(weto(i,j-1,k)*(THO(i,j-1,k)-THO(i,j,k))/dlyv(i,j-1)             &
     &+weto(i,j+1,k)*(THO(i,j+1,k)-THO(i,j,k))/dlyv(i,j))/dlyp(i,j))
!
      S1O(i,j,k)=weto(i,j,k)*(                                          &
     & (weto(i-1,j,k)*(SAO(i-1,j,k)-SAO(i,j,k))/dlxu(i-1,j)             &
     &+weto(i+1,j,k)*(SAO(i+1,j,k)-SAO(i,j,k))/dlxu(i,j))/dlxp(i,j)     &
     &+(weto(i,j-1,k)*(SAO(i,j-1,k)-SAO(i,j,k))/dlyv(i,j-1)             &
     &+weto(i,j+1,k)*(SAO(i,j+1,k)-SAO(i,j,k))/dlyv(i,j))/dlyp(i,j))
      enddo
      enddo
      enddo

!$OMP SINGLE
      CALL bounds_exch(T1O)
      CALL bounds_exch(S1O)
!$OMP END SINGLE

!$OMP DO
      do k=1,ke
      do j=2,je1
      do i=2,ie1
      THO(i,j,k)=THO(i,j,k)                                             &
     &- (weto(i,j,k)/((ddpo(i,j,k)+almzer                               &
     & +zzsurf(k)*(zo(i,j)-RHOICWA*sictho(i,j)-RHOSNWA*SICSNO(I,J)))    &
     &    *dlxp(i,j)*dlyp(i,j)))*dzw(k)*(                               &
     & weto(i-1,j,k)*aulx(i-1,j)*dlyu(i-1,j)*                           &
     &      (T1O(i-1,j,k)-T1O(i,j,k))/dlxu(i-1,j)                       &
     &+weto(i+1,j,k)*aulx(i,j)*dlyu(i,j)*                               &
     &       (T1O(i+1,j,k)-T1O(i,j,k))/dlxu(i,j)                        &
     &+weto(i,j-1,k)*auly(i,j-1)*dlxv(i,j-1)*                           &
     &       (T1O(i,j-1,k)-T1O(i,j,k))/dlyv(i,j-1)                      &
     &+weto(i,j+1,k)*auly(i,j)*dlxv(i,j)*                               &
     &       (T1O(i,j+1,k)-T1O(i,j,k))/dlyv(i,j))
      SAO(i,j,k)=SAO(i,j,k)                                             &
     &- (weto(i,j,k)/((ddpo(i,j,k)+almzer                               &
     & +zzsurf(k)*(zo(i,j)-RHOICWA*sictho(i,j)-RHOSNWA*SICSNO(I,J)))    &
     &    *dlxp(i,j)*dlyp(i,j)))*dzw(k)*(                               &
     & weto(i-1,j,k)*aulx(i-1,j)*dlyu(i-1,j)*                           &
     &      (S1O(i-1,j,k)-S1O(i,j,k))/dlxu(i-1,j)                       &
     &+weto(i+1,j,k)*aulx(i,j)*dlyu(i,j)*                               &
     &       (S1O(i+1,j,k)-S1O(i,j,k))/dlxu(i,j)                        &
     &+weto(i,j-1,k)*auly(i,j-1)*dlxv(i,j-1)*                           &
     &       (S1O(i,j-1,k)-S1O(i,j,k))/dlyv(i,j-1)                      &
     &+weto(i,j+1,k)*auly(i,j)*dlxv(i,j)*                               &
     &       (S1O(i,j+1,k)-S1O(i,j,k))/dlyv(i,j))
      enddo
      enddo
      enddo

!$OMP SINGLE
      CALL bounds_exch(THO)
      CALL bounds_exch(SAO)
!$OMP END SINGLE

      write(0,*)'done biharmonic: ',i,j,k
      ENDIF ! AULAPTS > ALZMER

       IF (I3DREST .GT. 0) THEN
       fakres=dt/(86400.*30.*4.)
!$OMP SINGLE
       write(io_stdout,*)'3D-RESTORING!!!!'
!$OMP END SINGLE
!$OMP DO
       do 2611 k=1,ke
       do 2611 j=1,je
       do 2611 i=1,ie
       if(ltlev(i,j,k))                                                 &
     &   tho(i,j,k)=tho(i,j,k)+fakres*(tlevi(i,j,k)-tho(i,j,k))
       if(lslev(i,j,k))                                                 &
     &   sao(i,j,k)=sao(i,j,k)+fakres*(slevi(i,j,k)-sao(i,j,k))
2611   continue
       ENDIF !I3DREST
!
!      schliesse imal
!
!$OMP END PARALLEL
      ENDDO ! IMAL
      RETURN
!
      END

