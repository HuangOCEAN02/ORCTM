      SUBROUTINE OCTDIFF_BASE
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_UNITS
      USE MO_OCTDIFF
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
      INTEGER i, j, k, imal
      REAL scale, slcut, stabmin, cutdel

      REAL AHX(IE,JE),AHY(IE,JE),AHXD(IE,JE),AHYD(IE,JE),ZZSURF(KE)
      REAL RHO_diff, XTERM,ZTERM,YTERM
!
      DO K=1,KE
       ZZSURF(K)=0.
      ENDDO
      ZZSURF(1)=1.
!
!     MINIMUM FOR STABILITY
!
      STABMIN=1.E-7
      CUTDEL=0.02
!
      DO IMAL=1,1

!$OMP PARALLEL PRIVATE(i,j,k,SCALE,SLCUT,RHO_diff,XTERM,ZTERM,YTERM)

      IF(AH00.GT.ALMZER)THEN

!$OMP DO
      DO J=1,JE
       DO I=1,IE
        AHX(I,J)=AH00*MAX(DLXU(I,J),DLYU(I,J))
        AHY(I,J)=AH00*MAX(DLYU(I,J),DLYV(I,J))
        AHXD(I,J)=AHX(I,J)*DT*DLYU(I,J)
        AHYD(I,J)=AHY(I,J)*DT*DLXV(I,J)
       ENDDO
      ENDDO

!$OMP DO
      DO K=1,KE                  
       DO J=1,JE 
        DO I=1,IE                
        vol_term(i,j,k)=DLXP(I,J)*DLYP(I,J)                             &
     &  *(ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)-RHOSNWA*SICSNO(I,J))   &
     &                +DDPO(I,J,K)+(1.-WETO(I,J,K)))
        ENDDO 
        KO(J,K)=MAX(K-1,1)
        KU(J,K)=MIN(K+1,KE)
        ZSURF(J,K)=FLOAT(K-KO(J,K))
        ZBOTT(J,K)=FLOAT(KU(J,K)-K)
       ENDDO
      ENDDO

!X DIRECTION- tracer independent values
!$OMP DO
      DO J=2,JE1
       DO I=1,IE1
        DO K=1,KE 
        IF(AMSUO(I,J,K).GT.0.)THEN
#ifdef ISOPYK
        slcut=CUTDEL*DZW(K)**2/(AHX(I,J)*DT)
        RHO_diff=(RHOO(I+1,J,K)-RHOO(I,J,K))/DLXU(I,J)
        XTERM=0.25*AHXD(I,J)*DDUO(I,J,K)
        ZTERM=0.25*DT*AHX(I,J)
!
!       TRIANGLE LEFT,UPW
!       
        SLOPLOx(I,J,K)= RHO_diff/MAX(STABIO(I,J,K),STABMIN)
        scale=1./MAX(1.,(SLOPLOx(I,J,K)/slcut)**2)
        XCOEFF_LO(I,J,K)=XTERM*ZSURF(J,K)*scale
        ZCOEFF_LOx(I,J,K)=ZTERM*ZSURF(J,K)*scale*SLOPLOx(I,J,K)
!       
!       TRIANGLE LEFT DOWN
!        
         SLOPLUx(I,J,K)= RHO_diff/MAX(STABIO(I,J,KU(J,K)),STABMIN)
         scale=1./MAX(1.,(SLOPLUx(I,J,K)/slcut)**2)
         XCOEFF_LU(I,J,K)=XTERM*ZBOTT(J,K)*WETO(I,J,KU(J,K))*scale
         ZCOEFF_LUx(I,J,K)=ZTERM*ZBOTT(J,K)*scale*SLOPLUx(I,J,K)*       &
     &                     WETO(I,J,KU(J,K))
!       
!       TRIANGLE RIGHT,UPW
!       
        SLOPROx(I,J,K)=RHO_diff/MAX(STABIO(I+1,J,K),STABMIN)
        scale=1./MAX(1.,(SLOPROx(I,J,K)/slcut)**2)
        XCOEFF_RO(I,J,K)=XTERM*ZSURF(J,K)*scale
        ZCOEFF_ROx(I,J,K)=ZTERM*ZSURF(J,K)*scale*SLOPROx(I,J,K)
!
!       TRIANGLE RIGHT DOWN
!
         SLOPRUx(I,J,K)=RHO_diff/MAX(STABIO(I+1,J,KU(J,K)),STABMIN)
         scale=1./MAX(1.,(SLOPRUx(I,J,K)/slcut)**2)
         XCOEFF_RU(I,J,K)=XTERM*ZBOTT(J,K)*WETO(I+1,J,KU(J,K))*scale
         ZCOEFF_RUx(I,J,K)=ZTERM*ZBOTT(J,K)*scale*SLOPRUx(I,J,K)*       &
     &                     WETO(I+1,J,KU(J,K))
#else
        XFLUX(I,J,K)=AHXD(I,J)*DDUO(I,J,K)
#endif
        ENDIF
        ENDDO
       ENDDO
      ENDDO
!Y DIRECTION - Tracer independent values
!$OMP DO
      DO J=1,JE1
       DO I=2,IE1
        DO K=1,KE
        IF(AMSUE(I,J,K).GT.0.)THEN
#ifdef ISOPYK
        SCALE=1.
        slcut=CUTDEL*DZW(K)**2/(AHY(I,J)*DT)
        RHO_diff=(RHOO(I,J+1,K)-RHOO(I,J,K))/DLYV(I,J)
        YTERM=0.25*AHYD(I,J)*DDUE(I,J,K)
        ZTERM=0.25*DT*AHY(I,J)
!
!       TRIANGLE LEFT,UPW
!
        SLOPLOy(I,J,K)= RHO_diff/MAX(STABIO(I,J,K),STABMIN)
        scale=1./MAX(1.,(SLOPLOy(I,J,K)/slcut)**2)
        YCOEFF_LO(I,J,K)=YTERM*ZSURF(J,K)*scale
        ZCOEFF_LOy(I,J,K)=ZTERM*ZSURF(J,K)*scale*SLOPLOy(I,J,K)
!
!       TRIANGLE LEFT DOWN
!
         SLOPLUy(I,J,K)= RHO_diff/MAX(STABIO(I,J,KU(J,K)),STABMIN)
         scale=1./MAX(1.,(SLOPLUy(I,J,K)/slcut)**2)
         YCOEFF_LU(I,J,K)=YTERM*ZBOTT(J,K)*WETO(I,J,KU(J,K))*scale
         ZCOEFF_LUy(I,J,K)=ZTERM*ZBOTT(J,K)*scale*SLOPLUy(I,J,K)*       &
     &                     WETO(I,J,KU(J,K))
!
!       TRIANGLE RIGHT,UPW
!
        SLOPROy(I,J,K)=RHO_diff/MAX(STABIO(I,J+1,K),STABMIN)
        scale=1./MAX(1.,(SLOPROy(I,J,K)/slcut)**2)
        YCOEFF_RO(I,J,K)=YTERM*ZSURF(J,K)*scale
        ZCOEFF_ROy(I,J,K)=ZTERM*ZSURF(J,K)*scale*SLOPROy(I,J,K)

!
!       TRIANGLE RIGHT DOWN
!
         SLOPRUy(I,J,K)=RHO_diff/MAX(STABIO(I,J+1,KU(J,K)),STABMIN)
         scale=1./MAX(1.,(SLOPRUy(I,J,K)/slcut)**2)
         YCOEFF_RU(I,J,K)=YTERM*ZBOTT(J,K)*WETO(I,J+1,KU(J,K))*scale
         ZCOEFF_RUy(I,J,K)=ZTERM*ZBOTT(J,K)*scale*SLOPRUy(I,J,K)*       &
     &                    WETO(I,J+1,KU(J,K))
#else
        YFLUX(I,J,K)=AHYD(I,J)*DDUE(I,J,K)
#endif
        ENDIF
        ENDDO
       ENDDO
      ENDDO
 
      ENDIF ! AH00.GT.ALMZER

! VERTICAL DIFFUSION - tracer independent values
!     IMPLICIT VERTICAL DIFFUSION

!$OMP DO
      DO J=2,JE1
      DO I=2,IE1
!
      DO K=1,KE
! INCLUDE ACTUAL LEVEL THICKNESS
      TRIDSY(I,J,K,1) = - DT*DVO(I,J,K)*WETO(I,J,K)*DI(K)               &
     & /(DDPO(I,J,K)+ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)             &
     &                      -RHOSNWA*SICSNO(I,J))+ALMZER)
      TRIDSY(I,J,K,3) = - DT*DVO(I,J,K+1) * DI(K+1)                     &
     & /(DDPO(I,J,K)+ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)             &
     &                      -RHOSNWA*SICSNO(I,J))+ALMZER)
      TRIDSY(I,J,K,2) = 1. - TRIDSY(I,J,K,1) - TRIDSY(I,J,K,3)
      END DO
!
      DO K=2,KE
      TRIDSY(I,J,K-1,1) = TRIDSY(I,J,K,1) / TRIDSY(I,J,K-1,2)
      TRIDSY(I,J,K,2)   = TRIDSY(I,J,K,2)                               &
     &  - TRIDSY(I,J,K-1,3) * TRIDSY(I,J,K,1) / TRIDSY(I,J,K-1,2)
      END DO
!
      END DO
      END DO
!
! BIHARMONIC DIFFUSION - tracer independent values
!
       IF(AULAPTS.GT.ALMZER)THEN
!$OMP DO
       do j=1,je
        do i=1,ie
#ifndef AULREDSC
        AULX(I,J)=AULAPTS*DLXU(I,J)**4
        AULY(I,J)=AULAPTS*DLYV(I,J)**4
#else
        AULX(I,J)=1.E4*AULAPTS*DLXU(I,J)**3
        AULY(I,J)=1.E4*AULAPTS*DLYV(I,J)**3
#endif
        enddo
       enddo
       ENDIF
!
!      schliesse imal
!
!$OMP END PARALLEL
      ENDDO ! IMAL
      RETURN
!
      END
