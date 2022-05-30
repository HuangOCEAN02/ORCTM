      SUBROUTINE GROWTH(ALAT,HO,AO,HSNO,HN,AN,HSNN,                     &
     &                  TICE,RPREC,TAIR,TD,ACL,PA,UG,SLN,               &
     &                  QS,QT,DZ,Z,DT,FW,SH,OM,                         &
     &                  QNSW,QNLW,QNSE,QNLA,PRECN,PRECH,HEATABS,SWSUM)
!=======================================================================
!
!     PROGRAMMED BY:
!     --------------
!     B.OWENS, P.LEMKE
!
!     MODIFIED BY:
!     ------------
!     A. STOESSEL                 MPI,  HAMBURG                     1989
!     S. LEGUTKE                  DKRZ, HAMBURG               24.12.1995
!
!        - ONLY SURFACE RESIDUAL HEAT USED FOR SNOW MELT
!
!        - SNOW --> ICE CONVERSION BECOMES SNOW --> WATER CONVERSION
!          IF ICE THICKNESS > HSNOTOICE 
!
!     S. LEGUTKE                  DKRZ, HAMBURG               24.06.1997
!        - ADJUST COMPACTNESS IN ORDER TO HAVE MINIMUM ICE THICKNESS SICTHMIN
!                
!        - COMPILE OPTION FOR FORCING WITH FLUXES INSTEAD OF SURFACE FIELDS 
!          (FLUXES, USED WITH THE COUPLED MODEL)
!
!     UWE MIKOLAJEWICZ     5/99
!        -INCLUDE EVAPORATION (FROM LATENT HEAT) IN FORING WITH SURFACE FIELDS
!
!     S. VENZKE                    MPI, HAMBURG               31.08.1999
!         -ADJUST COMPILE OPTION FOR FORCING WITH FLUXES INSTEAD OF SURFACE 
!          FIELDS
!
!     UWE                                                     28.06.00
!         INCLUDE SW-PENETRATION UNDER SEA ICE
!
!     S. LEGUTKE                  DKRZ, HAMBURG               04.08.2001
!        - At exit PRECN [m/s] contains all atmospheric freshwater fluxes
!          relevant for tracer concentration (i.e. snowmelt, river runoff,
!          P-E over open water, but no snow sublimation). PRECN also
!          contains snow->water conversion for cases of very thick ice.
!        - At exit BRINE [m/s] contains all changes of mean ice thickness
!          (eq. water depth) relevant for tracer concentration.
!        - BRINE(,) -> COMMON COMMOAU3 for use in marine bgc tracer dilution.
!
!
!     PURPOSE:
!     --------
!     UPDATE OF ICE THICKNESS, COMPACTNESS, SNOW DEPTH, UPPER OCEAN
!     TEMPERATURE AND SALINITY DUE TO ATMOSPHERIC HEAT AND FRESHWATER FLUX.
!
!     METHOD:
!     -------
!     GROWTH RATES OF THIN/THICK ICE AND SNOW DEPTH DUE TO ATMOSPHERIC 
!     FLUXES OF HEAT AND FRESHWATER ARE COMPUTED BASED ON THICKNESS
!     AND COMPACTNESS BEFORE THE UPDATE BY ADVECTIVE CHANGES IN SBR OCICE,
!     I.E. H(,,1) AND A(,,1) (I.E. HO AND AO) ARE USED IN BUDGET.
!     THESE GROWTH RATES ARE THEN USED TO UPDATE THICKNESS AND COMPACTNESS
!     AS OBTAINED IN SBR OCICE, I.E. H(,,2) AND A(,,2) (HN AND AN)
!     ARE UPDATED IN GROWTH.
!     NOTE: SURFACE ELEVATION Z CONTAINS ICE AND SNOW LAYER FOR DYNAMICS,
!     THUS ALL PRECIPITATION HAS TO BE ADDED TO Z, EVEN SNOWFALL.
!     FOR MIXING PROCESSES ICE+SNOW DRAFT HAS TO BE SUBTRACTED.
!
!SV 21.04.99 INCLUDED OPTION FOR FORCING WITH FLUXES

!C     INTERFACE:
!     ----------
!     INPUT:
!     *HO*       OLD ICE THICKNESS                      [M]
!     *AO*       OLD ICE COMPACTNESS                    [FRAC.]
!     *HSNO*     OLD SNOW DEPTH                         [M]
!     *RPREC*    ATMOSPHERIC P-E                        [M/S]
!     *TAIR*     AIR TEMPERATURES                       [DEG C]
!UWE     *SAIR*     SURFACE TEMPERATURES                   [K]
!     *TD*       DEW POINT TEMPERATURES                 [K]
!     *ACL*      FRACTIONAL CLOUD COVER                 [FRAC.]
!     *PA*       ATMOSPHERIC SURFACE PRESURE            [PA]
!     *UG*       WIND SPEED                             [M/S]
!UWE     *AOTX*     WIND STRESS, ZONAL                     [PA]
!UWE     *AOTY*     WIND STRESS, MERIDIONAL                [PA]
!     *SLN*      INCOMING SURFACE SOLAR RAD.            [W/M**2]
!     *DZ*       UNDISTURBED 1.LAYER THICK.             [M]
!     *DT*       TIME STEP                              [S]
!UWE     *AODFLX*   HEAT FLUX CORRECTION                   [W/M**2]
!     *OM*       LAND/SEA MASK                          [0/1]
!UWE     *AOHFLX*   ATMOSPHERIC HEAT FLUX                  [W/M**2]
!
!     OUTPUT:
!     *FW*       CHANGE OF ICE THICKNESS (MELT > 0)[M]
!                INCLUDING SNOW-TO-ICE CONVERSION
!     *SH*       ATMOS.hEAT FLUX (WITHOUT CORR.)   [W/M**2]
!     *PRECN*    P-E MODIFIED BY SNOW CHANGE       [M/S]
!     *QNSW*     ABSORBED SOLAR RADIATION          [W/M**2]
!     *QNLW*     OUTGOING LONGWAVE HEAT FLUX       [W/M**2]
!     *QNSE*     SENSIBLE HEAT FLUX                [W/M**2]
!     *QNLA*     LATENT HEAT FLUX                  [W/M**2]
!
!     CHANGED:
!     *Z*        SEA SURFACE ELEVATION                  [M]
!     *HN*       NEW ICE THICKNESS           [M]
!     *AN*       NEW ICE COMPACTNESS         [FRAC.]
!     *HSNN*     NEW SNOW DEPTH              [M]
!     *TICE*     ICE/SNOW SURFACE TEMPERATURE[DEG C]
!     *QS*       OCEANIC UPPER LAYER SALINITY[PSU]
!     *QT*       OCEANIC UPPER LAYER TEMP.   [DEG C]

!SVX 31.08.99
!
!     EXTERNALS:
!     ----------
!     OBUDGET: CALCULATES OPEN WATER ICE GROWTH RATES
!      BUDGET: CALCULATES ICE GROWTH RATES OVER ICE
!
!=======================================================================
!
!INCLUDE PARAM1
      USE MO_PARAM1
      USE MO_COMMOAU1 
      USE MO_COMMOAU3 
      USE MO_UNITS
!SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES

!
      DIMENSION HO(IE,JE),AO(IE,JE),HSNO(IE,JE)
      DIMENSION TAIR(IE,JE),RPREC(IE,JE)
      DIMENSION TD(IE,JE),ACL(IE,JE),PA(IE,JE)
!UWE      DIMENSION AOTX(IE,JE),AOTY(IE,JE),SAIR(IE,JE)
      DIMENSION SLN(IE,JE),Z(IE,JE)
      DIMENSION OM(IE,JE)
!UWE      DIMENSION AODFLX(IE,JE),AOHFLX(IE,JE)
!
      DIMENSION FW(IE,JE),SH(IE,JE),PRECN(IE,JE),PRECH(IE,JE)
      DIMENSION QNSW(IE,JE),QNLW(IE,JE),QNSE(IE,JE),QNLA(IE,JE)
!
      DIMENSION HN(IE,JE),AN(IE,JE),HSNN(IE,JE)
      DIMENSION TICE(IE,JE)
      DIMENSION QS(IE,JE),QT(IE,JE)
      DIMENSION ALAT(IE,JE)
!
!  LOCAL VARIABLES
!UWE      DIMENSION QTM(IE,JE) 
      DIMENSION TICM(IE,JE), TICA(IE,JE)
      DIMENSION UG(IE,JE),HEATABS(IE,JE)
      DIMENSION FO(IE,JE)
      DIMENSION QH(IE,JE),HS(IE,JE),HT(IE,JE)
      DIMENSION HDRAFT(IE,JE)
      DIMENSION RA(IE,JE),RHS(IE,JE),RHB(IE,JE),RH(IE,JE)
      DIMENSION RHSA(IE,JE),RHBA(IE,JE)
      DIMENSION QHST(IE,JE),SN(IE,JE),QFM(IE,JE)
      DIMENSION H(IE,JE,2),A(IE,JE,2),HSN(IE,JE,2)
      DIMENSION QSW(IE,JE),QLW(IE,JE),QSE(IE,JE),QLA(IE,JE)
      DIMENSION TMP(IE,JE),TMP2(IE,JE),TMP3(IE,JE),TMP4(IE,JE)

!SVX 31.08.99    
!
!UWE     ADDITIONAL LOCAL FIELD
!
      DIMENSION TMPUWE(IE,JE),TMPUW2(IE,JE),TMPUW3(IE,JE),TMPUW4(IE,JE)


      L = IE
      M = JE
!
!-----------------------------------------------------------------------
!  OLD AND NEW TIME LEVEL INDEX
!-----------------------------------------------------------------------
!
      LRHS=1
      LNEW=2
!
!-----------------------------------------------------------------------
!  DETERMINE ICE DRAFT AND TOTAL LAYER THICKNESS:
!-----------------------------------------------------------------------
!
      DO 1 J=1,M
      DO 1 I=1,L
         FO(I,J)=0.0  

         QSW(I,J)=0.0  
         QLW(I,J)=0.0  
         QSE(I,J)=0.0
         QLA(I,J)= 0.0 

         H  (I,J,LNEW) = HN  (I,J)
         A  (I,J,LNEW) = AN  (I,J)
         HSN(I,J,LNEW) = HSNN(I,J)
!
         FW (I,J)      = HN(I,J)
!
         H  (I,J,LRHS) = HO  (I,J)
         A  (I,J,LRHS) = AO  (I,J)
         HSN(I,J,LRHS) = HSNO(I,J)
!
         HDRAFT(I,J)=(RHOSNO*HSN(I,J,LNEW) + RHOICE*H(I,J,LNEW))/RHOWAT
!
         QH(I,J)= DZ + Z(I,J) - HDRAFT(I,J)
!
 1    CONTINUE
!
!-----------------------------------------------------------------------
!  NEW-ICE GROWTH RATE:
!-----------------------------------------------------------------------
!
!SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES
!
      CALL OBUDGET(ALAT,QT,TAIR,TD,ACL,PA,UG,SLN,FO                     &
     &            ,OM,QSW,QLW,QSE,QLA)
!
      VAPLI=1./VAPL

      DO 11 J=1,M
      DO 11 I=1,L
!
!SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES
         FO(I,J)=FO(I,J)+SWSUM*QSW(I,J)/CLB
         HEATABS(I,J)=SWSUM*QSW(I,J)*(1.-A(I,J,LRHS))

         QNSW(I,J) = QSW(I,J)*(1.-A(I,J,LRHS))
         QNLW(I,J) = QLW(I,J)*(1.-A(I,J,LRHS))
         QNSE(I,J) = QSE(I,J)*(1.-A(I,J,LRHS))
         QNLA(I,J) = QLA(I,J)*(1.-A(I,J,LRHS))
!UWE     STORE EVAPORATION IN TMPUW3
!
         TMPUW3(I,J)=DT*OM(I,J)*QLA(I,J)*(1.-A(I,J,LRHS))*VAPLI
!
!-----------------------------------------------------------------------
!  THICK ICE GROWTH RATES.
!-----------------------------------------------------------------------
!SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES
!-----------------------------------------------------------------------
!  CALCULATE EFFECTIVE ICE THICKNESS FOR CONDUCTION TERM IN BUDGET:
!               FIRST MAKE SURE WE HAVE NON-ZERO COMPACTNESS.
!               INCLUDE SNOW THICKNESS BECAUSE OF CONDUCTION EFFECT.
!               MAKE SURE TO HAVE NON-ZERO EFFECTIVE ICE THICKNESS.
!               INITIALISE MEAN ICE GROWTH RH AND TEMPERATURE (TICM).
!               STORE EFFECTIVE ICE THICKNESS IN ARRAY SN.
!-----------------------------------------------------------------------
!
         TMP (I,J) = ( H(I,J,LRHS) + HSN(I,J,LRHS)*CON/CONSN )          &
     &            /MAX(A(I,J,LRHS),ARMIN)
         RHS(I,J) = 0.0
         RHB(I,J) = 0.0
         RHSA(I,J) = 0.0
         RHBA(I,J) = 0.0
         TICM(I,J) = 0.0
         SN(I,J)=MAX(TMP(I,J),HMIN)
!
!-----------------------------------------------------------------------
!     SET ALBEDO (TMP2) ACCORDING TO PRESENCE OF SNOW  (TMP4)
!                                 TO MELTING CONDITIONS(TMP3).
!                                 USEMEAN PREVIOUS ICE TEMP..
!-----------------------------------------------------------------------
!
!UWE    INCLUDE FINITE VALUE FOR SNOW THICKNESS IN ALBEDO CALCULATION!
!
         TMP4(I,J) = (0.5-SIGN(0.5,1.E-2-HSN(I,J,LRHS)))
         TMP3(I,J) = (0.5+SIGN(0.5,TICE(I,J)))
         TMP2(I,J) = TMP4(I,J) *(TMP3(I,J)*ALBSNM+(1.-TMP3(I,J))*ALBSN) &
     &          +(1.-TMP4(I,J))*(TMP3(I,J)*ALBM  +(1.-TMP3(I,J))*ALBI )
!
!-----------------------------------------------------------------------
!  CALCULATE GROWTH RATES FOR THICK ICE.
!  DO ONCE FOR EACH ICE THICKNESS CATEGORIE ( ICELEV).
!  USE THE SAME MEAN ICE SURFACE TEMPERATURE OF THE PREVIOUS TIME
!  STEP FOR ALL THICKNESS CATEGORIES (LOOP 6).
!-----------------------------------------------------------------------
!
!UWE   INITIALIZE TEMPORARY FILEDS FOR SNOW EVAPORATION
!
        TMPUWE(I,J)=0.
        TMPUW2(I,J)=0.
        TMPUW4(I,J)=0.
 11   CONTINUE
!
!
      ICELEV=1
      ANZLEVI = 1./FLOAT(ICELEV)
!
      DO K = 1,ICELEV
!
      DO  J=1,M
      DO  I=1,L
!
         TMP (I,J) = (2*K-1)*SN(I,J)*ANZLEVI
         TICA(I,J) = TICE(I,J)
!
      ENDDO
      ENDDO
!
      CALL BUDGET(ALAT,RHSA,RHBA,TICA,LRHS,A,                           &
     &    TAIR,TD,ACL,PA,UG,SLN,OM,TMP,TMP2,TMP4,QSW,QLW,QSE,QLA)

!
!UWE   INCLUDE SUBLIMATION OF SNOW AND ICE FROM EVAPORATION
!
      SUBRSNI=1./(SUBL*RHOSNO)
      SUBRICI=1./(SUBL*RHOICE)
      RWRI=RHOSNO/RHOICE
      SUBLI=1./SUBL
!
      DO 13 J=1,M
      DO 13 I=1,L
         TICM(I,J) = TICM(I,J)+TICA(I,J)*ANZLEVI
         RHS (I,J) = RHS (I,J)+RHSA(I,J)*ANZLEVI
         RHB (I,J) = RHB (I,J)+RHBA(I,J)*ANZLEVI
         QNSW(I,J) = QNSW(I,J)+QSW(I,J)*A(I,J,LRHS)*ANZLEVI
         QNSE(I,J) = QNSE(I,J)+QSE(I,J)*A(I,J,LRHS)*ANZLEVI
         QNLW(I,J) = QNLW(I,J)+QLW(I,J)*A(I,J,LRHS)*ANZLEVI
         QNLA(I,J) = QNLA(I,J)+QLA(I,J)*A(I,J,LRHS)*ANZLEVI
!
!
!     TMPUWE CONTAINS AMOUNT OF SNOW AND ICE EVAPORATED DURING THIS TIME STEP
!     AVERAGED OVER GRID BOX   [KG/M**2]
!   
        TMPUWE(I,J)=TMPUWE(I,J)+A(I,J,LRHS)*QLA(I,J)                    &
     &        *ANZLEVI*OM(I,J)*SUBLI*DT
 13   CONTINUE
!
      ENDDO
!
!UWE        COMPUTE SUM OF SNOW AND ICE THICKNESS EVAPORATED AWAY.
!        IF NOT SUFFICIENT AVAILABLE, USE THE ENERGY TO EVAPORATE WATER
!        STORED IN TMPUW3
!
      VAPRHI=1./(VAPL*RHOWAT)
!
        DO J=1,M
         DO I=1,L
!SV 31.08.99 INCLUDED OPTION FOR FORCING WITH FLUXES
          TMPUW2(I,J)=-MIN(HSN(I,J,LNEW)*RHOSNO                         &
     &               +RHOICE*H(I,J,LNEW)                                &
     &            ,-TMPUWE(I,J))
!
          TMPUW3(I,J)=TMPUW3(I,J)+(TMPUWE(I,J)-TMPUW2(I,J))*SUBL/VAPL
          TMPUW4(I,J)=MIN(0.,TMPUW2(I,J)+HSN(I,J,LNEW)*RHOSNO)
         HSN(I,J,LNEW)=MAX(0.,HSN(I,J,LNEW)+TMPUWE(I,J)/RHOSNO)
         H(I,J,LNEW)=MAX(0.,H(I,J,LNEW)+TMPUW4(I,J)/RHOICE)
!
!       TMPUW2 CONTAINS SUBLIMATION OF ICE IN KG/M**2
!       TMPUW3 CONTAINS EVAPORATION OVER ICE FREE WATER IN KG/M**2
!
!-----------------------------------------------------------------------
!  STORE MEAN ICE TEMPERATURE ON TICE.
!  DETERMINE THERMODYNAMIC TOTAL ICE THICKNESS CHANGE (SH)
!  FOR CONTINUITY EQUATION.
!  MULTIPLY THICK ICE GROWTH RATE WITH COMPACTNESS A(,,LRHS) (SNOW DEPTH
!  IS MEAN OVER GRID CELL AREA).
!  CONVERT P-E OVER ICE INTO SNOWFALL IF TAIR < 0.
!  ADD SNOWFALL TO SNOW LAYER OF NEW TIME STEP.
!  SUBTRACT SNOWFALL FROM PRECIPITATION INTO MIXED LAYER.
!-----------------------------------------------------------------------
!
         RA (I,J)  = FO (I,J)*DT
         RHS(I,J)  = RHS(I,J)*DT
         RHB(I,J)  = RHB(I,J)*DT
         RH(I,J)  = RHS(I,J)+RHB(I,J)
!
!
         SH(I,J)       = RH(I,J)*A(I,J,LRHS) + RA(I,J)*(1.-A(I,J,LRHS))
!
         TICE(I,J)     = TICM(I,J)
         TMP3(I,J)= (0.5-SIGN(0.5,TAIR(I,J)))*A(I,J,LRHS)*FLOAT(ISNFLG)
         PRECN(I,J)    = ((1.-TMP3(I,J))* RPREC(I,J))*DT*OM(I,J)
         HSN(I,J,LNEW) = HSN(I,J,LNEW)                                  &
     &                  +    TMP3(I,J) * RPREC(I,J)*DT*OM(I,J)          &
     &                                             *RHOWAT/RHOSNO
         Z(I,J)        = Z(I,J) + RPREC(I,J)*DT*OM(I,J)
!
!UWE     ADD EVAPORATION TO SEA LEVEL
!
         Z(I,J)=Z(I,J)+OM(I,J)*(TMPUW2(I,J)+TMPUW3(I,J))/RHOWAT
         PRECN(I,J)=PRECN(I,J)+TMPUW3(I,J)/RHOWAT
!
         TMP(I,J)  = RHS(I,J)*A(I,J,LRHS)
!
!UWEX 10   CONTINUE
!
!
!-----------------------------------------------------------------------
!  AT MELTING CONDITIONS AT THE SURFACE (RHS<0), FIRST MELT SNOW 
!  MAKE SURE WE DO NOT END UP WITH NEGATIVE SNOW THICKNESS.
!-----------------------------------------------------------------------
       TMP2(I,J)=HSN(I,J,LNEW)+MIN(TMP(I,J),0.)*RHOICE/RHOSNO
!
!-----------------------------------------------------------------------
!  MODIFY ICE GROWTH TO ACCOUNT FOR HEAT LOSS BY SNOW MELT ( WITH CLO).
!  MODIFY PRECIPITATION TO ACCOUNT FOR SNOW MELT.
!  COMPUTE ICE MASS + HEAT STORAGE + ATMOS.hEAT FLUX (=-TOTAL HEAT) IN 
!  UPPER LAYER.
!  NOTE: IT IS ASSUMED THAT ALL PRECIPITATED WATER (P-E AND SNOW MELT)
!  HAS 2M-AIR TEMPERATURE (COMPILE OPTION FLUXES ONLY, SINCE IN THE
!  COUPLED VERSION WE DO NOT WANT TO TRANSFER 2M TEMPERATURE).
!-----------------------------------------------------------------------
       SN(I,J)=MAX(TMP2(I,J),0.)
         TMP2(I,J)    = HSN(I,J,LNEW)-SN(I,J)
         HSN(I,J,LNEW)= SN(I,J)
!
         RH(I,J)      = SH(I,J)   +TMP2(I,J)*RHOSNO/RHOICE
!
         PRECN(I,J)   = PRECN(I,J)+TMP2(I,J)*RHOSNO/RHOWAT
!
         QHST(I,J)    =  H(I,J,LNEW)                                    &
     &                 -( (QT  (I,J)-TFREZ)*QH(I,J)                     &
     &                   +(TAIR(I,J)-TFREZ)*PRECN(I,J)                  &
     &                   ) * CC/CLB*OM(I,J)                             &
     &                  + RH(I,J)
!
!-----------------------------------------------------------------------
!  WHEN NO ICE IS LEFT (QHST<0) MELT SNOW FIRST.
!-----------------------------------------------------------------------
        TMP(I,J)=MAX(QHST(I,J),0.)
        TMP3(I,J)=MIN(QHST(I,J),0.)
         SN(I,J)=SN(I,J)+TMP3(I,J)*RHOICE/RHOSNO
!-----------------------------------------------------------------------
!  UPDATE FRESH WATER FLUX INTO UPPER OCEAN LAYER (PRECN) WITH 
!         ADDITIONAL SNOW MELT.
!  UPDATE HEAT STORAGE IN UPPER OCEAN (QHST) DUE TO ADDITIONAL SNOW MELT.
!
!  UPDATE SNOW DEPTH WITH ADDITIONAL SNOW MELT (NOTE: IF SNOW REMAINS
!         WITHOUT ICE BENEATH, IT WILL BE TRANSFORMED INTO ICE (LOOP 141).
!UWE REMOVED C  QTM IS HEAT CONTENT IN UPPER OCEAN LAYER.
! 
!  H-QFM IS ICE THICKNESS CHANGE: > 0 ==> MELTING; < 0 ==> FREEZING.
!
!  UPDATE OCEANIC LAYER THICKNESS DUE TO ADDITIONAL SNOW MELT.
!-----------------------------------------------------------------------
       TMP2(I,J)=MIN(SN(I,J),0.)
       SN(I,J)=MAX(SN(I,J),0.)
       PRECN(I,J)= PRECN(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOSNO/RHOWAT
!
       QHST(I,J) = QHST(I,J)+(HSN(I,J,LNEW)-SN(I,J))*RHOSNO/RHOICE      &
     &              - (HSN(I,J,LNEW)-SN(I,J))*RHOSNO/RHOWAT             &
!UWE USE 0. INSTEAD OF TAIR FOR SNOW MELT
     &                                  *(0.-TFREZ)*CC/CLB

         QFM(I,J)      = H(I,J,LNEW)
         HSN(I,J,LNEW) = SN(I,J)
         RH(I,J)       =-RH(I,J)
!
!-----------------------------------------------------------------------
!  UPDATE ICE THICKNESS (EQ.9 IN OWENS & LEMKE 90)
!-----------------------------------------------------------------------
!
!  UPDATE SALT CONTENT AND
!         HEAT CONTENT OF UPPER OCEANIC LAYER.
!-----------------------------------------------------------------------
         H(I,J,LNEW)=MAX(QHST(I,J),0.)
         BRINE(I,J) = (QFM(I,J)-H(I,J,LNEW))*OM(I,J)
!
         HS(I,J)  = QS(I,J)*QH(I,J)                                     &
     &              +BRINE(I,J)*RHOICWA*SICE
!
         HT(I,J)  = -TMP2(I,J)*RHOSNO*CLB/(CC*RHOICE)
!
!-----------------------------------------------------------------------
!  UPDATE POTENTIAL TEMPERATURE AND SALINITY OF UPPER OCEANIC LAYER
!         AND LAYER THICKNESS FOR THERMODYNAMICS.
!-----------------------------------------------------------------------
         HDRAFT(I,J) =   RHOSNWA*HSN(I,J,LNEW)                          &
     &                  +RHOICWA*H(I,J,LNEW) 
!
         QH(I,J) = DZ + Z(I,J) - HDRAFT(I,J)
!
         QT(I,J) = (HT(I,J)/QH(I,J) + TFREZ)*    OM(I,J)                &
     &            + QT(I,J)                 *(1.-OM(I,J))
!
         QS(I,J) = HS(I,J)/QH(I,J)
!
!-----------------------------------------------------------------------
!  MAKE SURE WE DON'T TRY TO MELT MORE ICE THAN IS AVAILABLE:
!  NOTE: POSITIVE RH MEANS MELTING HERE.
!-----------------------------------------------------------------------
!
       RH(I,J)=-MIN(RH(I,J),H(I,J,LNEW))
       TMP3(I,J)=MAX(H(I,J,LRHS),HMIN)
       TMP2(I,J)=MIN(RH(I,J),0.)
       TMP(I,J)=MAX(RA(I,J),0.)
!-----------------------------------------------------------------------
!  UPDATE ICE COMPACTNESS (EQ.16 IN HIBLER 79)
!  MAKE SURE WE DO NOT DIVIDE BY 0 
!  IF MELTING THICK ICE, THEN EVALUATE THE MELTING TERM: TMP2
!  IF FREEZING THIN ICE, THEN EVALUATE THE FREEZING TERM: TMP.
!-----------------------------------------------------------------------
         RA(I,J)     = 0.5*TMP2(I,J) * A(I,J,LRHS)/TMP3(I,J)            &
! MODIFY LEADCLOSING IN CASE OF FREEZING
     &                +    TMP (I,J) *(1.-A(I,J,LRHS))*5.               &
#ifndef LEADCLOSE
     &      /(H0+4.*MAX(H0,H(I,J,LRHS)/MAX(A(I,J,LRHS),ARMIN)))
#else
     &      /(H0+3.*MAX(H0,H(I,J,LRHS)/MAX(A(I,J,LRHS),ARMIN)))
#endif
!
!
         A(I,J,LNEW) = A(I,J,LNEW) + RA(I,J)
!-----------------------------------------------------------------------
!  ENSURE THAT COMPACTNESS > 0 WHERE THERE IS ICE.
!  SET COMPACTNESS TO 0 WHERE THERE IS NO ICE.
!  COMPACTNESS IS NOT ALLOWED TO BECOME LARGER THAN ARMAX.
!  COMPACTNESS IS NOT ALLOWED TO BECOME LESS THAN 0.
!-----------------------------------------------------------------------
         IF( H(I,J,LNEW) .GT. 0. .AND. A(I,J,LNEW) .LE. 0.) THEN
            A(I,J,LNEW) = H(I,J,LNEW)/SICTHMIN
         ENDIF
!
         A(I,J,LNEW) = ( A(I,J,LNEW)*(0.5+SIGN(0.5, A(I,J,LNEW)      )) &
     &                              *(0.5-SIGN(0.5, A(I,J,LNEW)-ARMAX)) &
     &                  +      ARMAX*(0.5+SIGN(0.5, A(I,J,LNEW)-ARMAX)))&
     &                *(0.5-SIGN(0.5,-H(I,J,LNEW)))
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!  SNOW TO ICE CONVERSION:
!-----------------------------------------------------------------------
      IF ( ISNFLG .EQ. 1 ) THEN
!
!-----------------------------------------------------------------------
!     HSNTOICE IS  MAXIMUM ICE THICKNESS FOR WHICH SNOW --> ICE 
!     CONVERSION WILL BE DONE.
!-----------------------------------------------------------------------
!
         DO 140 J = 1,M
         DO 140 I = 1,L
!
            TMP3  (I,J)=(0.5-SIGN(0.5,H(I,J,LNEW)-HSNTOICE))
            TMP2(I,J)=MIN(HDRAFT(I,J),H(I,J,LNEW))
            TMP(I,J)=MAX(HDRAFT(I,J),H(I,J,LNEW))
 140     CONTINUE
!-----------------------------------------------------------------------
!  SNOW TO ICE CONVERSION:
!            IN CASE THE ICE SURFACE LIES IN A DEPTH DELTH BELOW THE 
!            WATER LINE BECAUSE OF THE SNOW LOAD , THE ICE THICKNESS 
!            IS INCREASED BY THIS DEPTH, AND AN EQUIVALENT AMOUNT (IN
!            IN HEAT) OF SNOW IS MELTED.
!            IF THERE IS NET SNOW ACCUMULATION AND NOT ENOUGH ICE MELT,
!            THIS CAN LEAD TO UNLIMITED ICE GROWTH. THUS IN CASE THE
!            ICE BECOMES TOO THICK LOCALLY, THE SNOW TO ICE CONVERSION
!            IS TRANSFORMED INTO A SNOW TO FRESHWATER CONVERSION THERE.
!  CLOSE THE SALINITY BALANCE.
!            NOTE:IF ICE IS FORMED, THE HEAT IS TAKEN FROM SNOW MELT,
!            I.E. NO IMPACT ON OCEAN TEMPERATURE; IF NO ICE IS FORMED
!            THE HEAT BALANCE IS NOT CLOSED.
!-----------------------------------------------------------------------
         DTI=1./DT
         DO 141 J=1,M
         DO 141 I=1,L
!
            DELTH         = HDRAFT(I,J)-TMP2(I,J)
!
            HSN(I,J,LNEW) = HSN(I,J,LNEW)-DELTH*RHOICE/RHOSNO
!
            H  (I,J,LNEW) =       TMP3(I,J) *TMP(I,J)                   &
     &                      +(1.0-TMP3(I,J))*H(I,J,LNEW)
!
            HDRAFT(I,J)   = ( RHOSNO*HSN(I,J,LNEW)                      &
     &                     +RHOICE*H(I,J,LNEW)  )/RHOWAT
!
            QH(I,J)       = DZ + Z(I,J) - HDRAFT(I,J)
!
            HS(I,J)       = HS(I,J)                                     &
     &                     -DELTH*SICE*RHOICE/RHOWAT*TMP3(I,J)
!
            QS(I,J)       = HS(I,J)/QH(I,J)
!
!UWE  CONVERT PRECN TO M/S
!
!CHANGES ACC. TO S.L.
         PRECH(I,J)=RPREC(I,J)+(TMPUW2(I,J)+TMPUW3(I,J))/(RHOWAT*DT)
            PRECN(I,J) = (PRECN(I,J)+DELTH*RHOICE/RHOWAT)/DT

            BRINE(I,J) = (BRINE(I,J)-DELTH*TMP3(I,J))*RHOICE/RHOWAT/DT
!
 141     CONTINUE
!
      ENDIF
!
!-----------------------------------------------------------------------
!  STORE HEAT FLUX (W/M**2) WITHOUT HEAT FLUX CORRECTION TERM IN SH.
!  FW IS TOTAL THERMODYNAMIC CHANGE OF ICE THICKNESS (M OF WATER).
!  ADJUST COMPACTNESS TO HAVE MINIMUM ICE THICKNESS SICTHMIN.
!-----------------------------------------------------------------------
!
      DO 142 J=1,M
      DO 142 I=1,L
!

         SH  (I,J) = -SH(I,J)*CLB/DT

         HN  (I,J) = H  (I,J,LNEW)
!
         AN  (I,J) = A(I,J,LNEW)                                        &
     &               *(0.5+SIGN(0.5,H(I,J,LNEW)/SICTHMIN-A(I,J,LNEW)))  &
     &              +H(I,J,LNEW)/SICTHMIN                               &
     &               *(0.5-SIGN(0.5,H(I,J,LNEW)/SICTHMIN-A(I,J,LNEW)))
!
         HSNN(I,J) = HSN(I,J,LNEW)
!
         FW  (I,J) =( FW(I,J)-H(I,J,LNEW) )*RHOICE/RHOWAT
!
  142 CONTINUE
!
      RETURN
      END
