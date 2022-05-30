      SUBROUTINE VARDRAG(DRAGL,DRAGS,K,SST,AIRT,HUMIDA,HUMIDO,UVABSOL)
#ifdef DASILVA
!===============================================================================
!
!     PURPOSE:  EVALUATION OF DRAG COEFFICIENT FOR MOMENTUM,SENSIBLE AND
!               LATENT HEAT
!
!     METHOD:   LARGE AND POND  (1981,1982)
!
!     AUTHOR:   JOSEF MAXIMILIAN OBERHUBER
!
!     DATE:     20. JULY 1991
!
!     MODIFIED FOR DA-SILVA-INTERCOMPARISON: 19. JUNE 1996 BY JMO
!
!     MODIFIED FOR INCLUDING IN HOPE:       02. JULY 1996 BY FR
!                                               (FRANK ROESKE)
!
!
!===============================================================================
!
      USE MO_PARAM1
      PARAMETER (RHALF=0.5)
      DIMENSION AIRT(IE*JE)  , UVABSOL(IE*JE) , UVARIAN(IE*JE)
      DIMENSION SST(IE*JE)   , DRAG(IE*JE)
      DIMENSION HUMIDO(IE*JE), HUMIDA(IE*JE)
      DIMENSION DRAGL(IE*JE) , DRAGS(IE*JE)
!
!---->DEFINE LOCAL VARIABLES
!
      DIMENSION HELP(IE*JE),PSIM(IE*JE),PSILS(IE*JE)
!
!---->INITIALIZE PARAMETERS
!
      SAVE CONSTS,CONSTL,REFHGT,CHARNCK,RKARMAN,GRAV
!
!---->INITIALIZE PARAMETERS
!
      DATA CONSTS  /.0327/   , CONSTL  /.0346/      , REFHGT /10./
!      DATA CHARNCK /.0064/    , RKARMAN /.4/         , GRAV   /9.806/
!
!---->TO CLOSE THE FRESH WATER BUDGET:
!
      DATA CHARNCK /.0187/    , RKARMAN /.4/         , GRAV   /9.806/
!
!---->DETERMINE STABILITY PARAMETER
!
!CC      WRITE(6,*) 'VARDRAG: K=',K
!
      DO I=1,K
        UVARIAN(I)=0.2*UVABSOL(I)
        DELTAV=UVABSOL(I)*(UVABSOL(I)+2.*UVARIAN(I))
        DELTAT=SST(I)-AIRT(I)
!       WRITE(6,*) 'I: ',I,' DELTAV: ',DELTAV,' UVABSOL: ',UVABSOL(I),
!     2           ' UVARIAN: ',UVARIAN(I),' SST: ',SST(I),
!     3           ' AIRT: ',AIRT(I)
!FR     ESA=611.*10**(7.5*(AIRT(I,J)-273.16)/(AIRT(I,J)-35.86))
!FR     ESO=611.*10**(7.5*(SST(I,J)-273.16)/(SST(I,J)-35.86))
!FR     HUMIDA=.662*ESA/(SLPRESS-(1.-.622)*ESA)*HUMID(I,J)
!FR     HUMIDO=.662*ESO/(SLPRESS-(1.-.622)*ESO)
        DELTAQ=HUMIDO(I)-HUMIDA(I)
        HELP(I)=GRAV*REFHGT/(CHARNCK*DELTAV)
!       WRITE(6,*) 'HELP: ',HELP(I)
        T0=AIRT(I)*(1.+1.7E-6*AIRT(I)*HUMIDA(I))
!       WRITE(6,*) 'T0: ',T0
        FLAG=(RHALF+SIGN(RHALF,DELTAT))
        ZOVERL=-70.*REFHGT*(DELTAT+2.5E-6*T0**2*DELTAQ)*(1.-FLAG)       &
     &             -100.*REFHGT*(DELTAT+1.7E-6*T0**2*DELTAQ)*FLAG
        ZOVERL=ZOVERL/(T0*DELTAV)
        X=SQRT(SQRT(1.+16.*ABS(ZOVERL)))
        PSIM(I)=-(1.-FLAG)*7.*ZOVERL                                    &
     &            +FLAG*(2.*LOG((1.+X)*.5)+LOG((1.+X**2)*.5)            &
     &                  +2.*(ATAN(1.)-ATAN(X)))
        PSILS(I)=(1.-FLAG)*PSIM(I)+FLAG*2.*LOG((1.+X**2)*.5)
!       WRITE(6,*) 'PSILS: ',PSILS(I)
        DRAG(I)=1.4E-3
      ENDDO
!
!---->SOLVE FOR DRAG COEFFICIENT
!
      SUM1=0.
      SUM2=0.
      DO ITER=1,5
       SUM=0.
!
       DO I=1,K
!
!      WRITE(6,*) 'DRAG: ',DRAG(I)
         XH = LOG(HELP(I)/DRAG(I))-PSIM(I)
         FVONX=SQRT(DRAG(I))*XH-RKARMAN
!      WRITE(6,*) 'XH-2/HELP: ',XH-2./HELP(I)
         DFDX=(2.*SQRT(DRAG(I)))/(XH-2./HELP(I))
         CM=DRAG(I)-FVONX*DFDX
         CM=MIN(MAX(CM,5.E-4),3.E-3)
         SUM=SUM+ABS(CM-DRAG(I))
         DRAG(I)=CM
!
!        XH = LOG(HELP(I)/DRAG(I))-PSIM(I)
!        FVONX=SQRT(DRAG(I))*XH-RKARMAN
!        DFDX=(2.*SQRT(DRAG(I)))/(XH-2.)
!        DRAG(I)=DRAG(I)-FVONX*DFDX
!        DRAG(I)=MAX(DRAG(I),5.E-4)
!        DRAG(I)=MIN(DRAG(I),3.E-3)
!        SUM=SUM+ABS(FVONX*DFDX)
!
       ENDDO
       SUM=SUM/FLOAT(K)
! RJ: We do not want that the behavior in a parallel run is different
!     because of different iteration numbers, we alway take
!     the maximum number of iterations:
!!!       IF(SUM.LT.1.E-6) GOTO 7
      ENDDO
    7 CONTINUE
!
!---->DETERMINE TRANSFER COEFFICIENTS FOR LATENT HEAT
!                                 AND FOR SENSIBLE HEAT
!
      DO I=1,K
         CMN=(RKARMAN/LOG(HELP(I)/DRAG(I)))**2
         CLN=CONSTL*RKARMAN/LOG(HELP(I)/DRAG(I))
         DRAGL(I)=CLN*SQRT(DRAG(I)/CMN)                                 &
     &            /(1.-CLN*PSILS(I)/(SQRT(CMN)*RKARMAN))
         CSN=CONSTS*RKARMAN/LOG(HELP(I)/DRAG(I))
         DRAGS(I)=CSN*SQRT(DRAG(I)/CMN)                                 &
     &            /(1.-CSN*PSILS(I)/(SQRT(CMN)*RKARMAN))
      ENDDO
      DO I=1,K
        DRAGL(I)=MAX(DRAGL(I),5.E-4)
        DRAGL(I)=MIN(DRAGL(I),3.E-3)
        DRAGS(I)=MAX(DRAGS(I),5.E-4)
        DRAGS(I)=MIN(DRAGS(I),3.E-3)
        SUM1=SUM1+DRAGS(I)
        SUM2=SUM2+DRAGL(I)
      ENDDO
!      WRITE(6,*)'VARDRAG: ',K,SUM1/FLOAT(K),SUM2/FLOAT(K)
#endif /*DASILVA*/
      RETURN
      END
