      SUBROUTINE BELEG_ZERO
!
!     BELEG_ZERO : INITIALIZES ALL FIELDS DEFINED in MO_COMMO1 by ZERO  
!---------------------------------------------------------------
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU2
      USE MO_COMMOAU3
      USE MO_ELICOM
      USE MO_PARA2
!
#ifdef SLOPECON_ADPO
      USE MO_COMMOBBL
#endif
!
#ifdef MEAN         
      USE MO_MEAN
#endif
!
!:: COMMON BLOCK GRID     
!:: THREE- DIM FIELDS
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         DDPO(I,J,K)=ZERO
         DPIO(I,J,K)=ZERO
         DDPSIO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=I_start,IE
         DDUO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=J_start,JE
        DO I=1,IE
         DDUE(I,J,K)=ZERO
        ENDDO
       ENDDO
      ENDDO
!:: TWO- DIM FIELDS
      DO J=1,JE
       DO I=1,IE
        DLXP(I,J)=ZERO
        DLYP(I,J)=ZERO
        DTDYO(I,J)=ZERO
        DTDXPO(I,J)=ZERO
        DEPTO(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=J_start,JE
       DO I=1,IE
        DLYV(I,J)=ZERO
        DLXV(I,J)=ZERO
        DTDXUE(I,J)=ZERO
        DPYE(I,J)=ZERO
        DEUTE(I,J)=ZERO
        DEUTIE(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=1,JE
       DO I=I_start,IE
        DLYU(I,J)=ZERO
        DLXU(I,J)=ZERO
        DTDXUO(I,J)=ZERO
        DPYO(I,J)=ZERO
        DEUTO(I,J)=ZERO
        DEUTIO(I,J)=ZERO
       ENDDO
      ENDDO
!:: ONE- DIM FIELDS
      DO K=1,KE     
       DWI(K)=ZERO
       DZW(K)=ZERO
      ENDDO
      DO K=1,KEP    
       DZ(K)=ZERO
       DI(K)=ZERO
       TIESTU(K)=ZERO
       TIESTW(K)=ZERO
      ENDDO
      DO J=1,JE
       SIFTFE(J)=ZERO
       SIFTFO(J)=ZERO
       COFTFE(J)=ZERO
       COFTFO(J)=ZERO
      ENDDO
!::
!:: COMMON BLOCK GEO
!:: TWO-DIMENSIONAL FIELDS
      DO J=1,JE
      DO I=1,IE
       ALAT(I,J)=ZERO
      ENDDO
      ENDDO
      DO J=1,JTO
      DO I=1,ITO
       GILA(I,J)=ZERO
       GIPH(I,J)=ZERO
      ENDDO
      ENDDO
!:: COMMON BLOCK MASK
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         WETO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=I_start,IE
         AMSUO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=J_start,JE
        DO I=1,IE
         AMSUE(I,J,K)=ZERO
        ENDDO
       ENDDO
      ENDDO
!:: COMMON BLOCK FIELDS3D
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         RHOO(I,J,K)=ZERO
         THO(I,J,K)=ZERO
         SAO(I,J,K)=ZERO
         PO(I,J,K)=ZERO
         T1O(I,J,K)=ZERO
         S1O(I,J,K)=ZERO
         STABIO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=I_start,IE
         UKO(I,J,K)=ZERO
         UK1O(I,J,K)=ZERO
         UOO(I,J,K)=ZERO
        ENDDO
       ENDDO
       DO J=J_start,JE
        DO I=1,IE
         VKE(I,J,K)=ZERO
         VK1E(I,J,K)=ZERO
         VOE(I,J,K)=ZERO
        ENDDO
       ENDDO
      ENDDO

      DO K=1,KEP
       DO J=1,JE
        DO I=1,IE
         DVO(I,J,K)=ZERO
         AVO(I,J,K)=ZERO
         WO(I,J,K)=ZERO
        ENDDO
       ENDDO
      ENDDO
!:: COMMON BLOCK FIELDS2D
      DO J=1,JE
       DO I=1,IE
        ZO(I,J)=ZERO
        Z1O(I,J)=ZERO
        SHELP(I,J)=ZERO
        THELP(I,J)=ZERO
        RHELP(I,J)=ZERO
        TLOW(I,J)=ZERO
        SLOW(I,J)=ZERO
        TSUP(I,J)=ZERO
        SSUP(I,J)=ZERO
        B1O(I,J)=ZERO
        NUM(I,J)=0   
        BRINE(I,J)=ZERO
        KCONDEP(I,J)=0
       ENDDO
      ENDDO
      DO J=1,JE
       DO I=I_start,IE
        U1O(I,J)=ZERO
        USO(I,J)=ZERO
        UCOS(I,J)=ZERO
        PXOIN(I,J)=ZERO
        UZO(I,J)=ZERO
        FTWOU(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=J_start,JE
       DO I=1,IE
        V1E(I,J)=ZERO
        VSE(I,J)=ZERO
        VCOS(I,J)=ZERO
        PYEIN(I,J)=ZERO
        VZE(I,J)=ZERO
        FTWOV(I,J)=ZERO
       ENDDO
      ENDDO
!:: COMMON BLOCK FIELDS1D
      DO K=1,KE
      DOTRT(K)=ZERO
      UPTRT(K)=ZERO
      PREFF(K)=ZERO
      ROREF(K)=ZERO
      ENDDO
      DO K=1,KEP
      WINTUR(K)=ZERO
      DBACKV(K)=ZERO
      ABACKV(K)=ZERO
      ENDDO
      DO i=1,ill
      SKAL(I)=ZERO
      ENDDO
      DO i=1,ilt
      B(I)=ZERO
      X(I)=ZERO
      ENDDO
!:: COMMON BLOCK HOPEICF
      DO J=1,JE
       DO I=1,IE
        SICTHO(I,J)=ZERO
        SICOMO(I,J)=ZERO
        SICOMU(I,J)=ZERO
        SICOMV(I,J)=ZERO
        SICOMP(I,J)=ZERO
        SICUO(I,J)=ZERO
        SICVE(I,J)=ZERO
        HIBDELO(I,J)=ZERO
        HIBZETO(I,J)=ZERO
        HIBETO(I,J)=ZERO
        HIBDELE(I,J)=ZERO
        HIBZETE(I,J)=ZERO
        HIBETE(I,J)=ZERO
        TAUWATU(I,J)=ZERO
        TAUWATV(I,J)=ZERO
        SICSNO(I,J)=ZERO
        FWO(I,J)=ZERO
        SHO(I,J)=ZERO
       ENDDO
      ENDDO
!::
      DO J=1,JE
       DO I=1,IE
        KBOT(I,J)=0   
        TAFO(I,J)=ZERO
        EMINPO(I,J)=ZERO
        FCLOU(I,J)=ZERO
        FSWR(I,J)=ZERO
        FU10(I,J)=ZERO
        FPREC(I,J)=ZERO
        FTDEW(I,J)=ZERO
        RELSAO(I,J)=ZERO
        RELTHO(I,J)=ZERO
        GIRIV(I,J)=ZERO
        RIVRUN(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=1,JE
       DO I=I_start,IE
        TXO(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=J_start,JE
       DO I=1,IE
        TYE(I,J)=ZERO
       ENDDO
      ENDDO
!
#ifdef GMBOLUS
       DO K=1,KEP
        DO J=1,JE
         DO I=1,IE
          WGO(I,J,K)=ZERO
         ENDDO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=1,IE
         BOLX(I,J)=ZERO
         BOLY(I,J)=ZERO
        ENDDO
       ENDDO
#endif /*GMBOLUS*/
!:: USE FILE MO_COMMO2
!:: COMMON BLOCK FORUWE2
      DO J=1,JE
       DO I=1,IE
        TAFO1(I,J)=ZERO
        FTDEW1(I,J)=ZERO
        FSWR1(I,J)=ZERO
        FSLP1(I,J)=ZERO
        FU101(I,J)=ZERO
        FCLOU1(I,J)=ZERO
        FPREC1(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=1,JE
       DO I=I_start,IE
        TXO1(I,J)=ZERO
       ENDDO
      ENDDO
      DO J=J_start,JE
       DO I=1,IE
        TYE1(I,J)=ZERO
       ENDDO
      ENDDO
!:: USE FILE MO_COMMOAU2
!:: COMMON bLOCK SIOBOD
       DO J=1,JE
        DO I=1,IE
         TICEO(I,J)=ZERO
         ACLO(I,J)=ZERO
         PAO(I,J)=ZERO
         RPRECO(I,J)=ZERO
         FRSE(I,J)=ZERO
         PRECO(I,J)=ZERO
         PRECH(I,J)=ZERO
         QSWO(I,J)=ZERO
         QLWO(I,J)=ZERO
         QSEO(I,J)=ZERO
         QLAO(I,J)=ZERO
         TAIRO(I,J)=ZERO
         TDO(I,J)=ZERO
         SICUDO(I,J)=ZERO
         SICVDE(I,J)=ZERO
        ENDDO
       ENDDO
!:: INCLUDE FILE PARA2.h
!:: COMMON BLOCK ITHILF
       DO J=1,JE
        DO I=1,IE
         FF(I,J)=ZERO
         XX(I,J)=ZERO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=I_start,IE
         UF(I,J)=ZERO
        ENDDO
       ENDDO
       DO J=J_start,JE
        DO I=1,IE
         VF(I,J)=ZERO
        ENDDO
       ENDDO
!:: INCLUDE FILE MO_ELICOM.h
!:: COMMON BLOCK DUMDUM
       DO J=1,ILL
        DO I =1,IMM
         PGL(I,J)=ZERO
        ENDDO
       ENDDO
#ifdef SLOPECON_ADPO
!:: INCLUDE FILE MO_COMMOBBL
!:: COMMON BLOCL SLOPE
       DO J=1,JE
        DO I=1,IE
         ALPX(I,J)=ZERO
         ALPY(I,J)=ZERO
         BBLFLX(I,J)=ZERO
         UBBL(I,J)=ZERO
         VBBL(I,J)=ZERO
         KUPUBBL(I,J)=0   
         KDWUBBL(I,J)=0   
         KUPVBBL(I,J)=0   
         KDWVBBL(I,J)=0   
        ENDDO
       ENDDO
#endif /*SLOPECON_ADPO*/
!::
#ifdef MEAN
!:: COMMON BLOCK MMEAN
       DO K=1,KE
        DO J=1,JE
         DO I=I_start,IE
          SUM_UKO(I,J,K)=ZERO
         ENDDO
        ENDDO
        DO J=J_start,JE
         DO I=1,IE
          SUM_VKE(I,J,K)=ZERO
         ENDDO
        ENDDO
        DO J=1,JE
         DO I=1,IE
         SUM_THO(I,J,K)=ZERO
         SUM_SAO(I,J,K)=ZERO
         SUM_WO(I,J,K)=ZERO
         SUM_PO(I,J,K)=ZERO
! PO added by peterspy
#ifdef DIFFDIAG
         SUM_AVO(I,J,K)=ZERO
         SUM_DVO(I,J,K)=ZERO
         SUM_WTMIX(I,J,K)=ZERO
         WTMIX(I,J,K)=ZERO
#endif /*DIFFDIAG*/

#ifdef NON_HYDROSTATIC
         SUM_PNH(I,J,K)=ZERO
         SUM_DIVG(I,J,K)=ZERO
#ifdef CG2D4NONHY             
         SUM_PNHBOT(I,J)=ZERO
         SUM_DIVGBOT(I,J)=ZERO
#endif /*CG2D4NONHY*/        
#endif /*NON_HYDROSTATIC*/

#ifdef GMBOLUS
         SUM_WGO(I,J,K)=ZERO
         SUM_BOLX(I,J)=ZERO
         SUM_BOLY(I,J)=ZERO
#endif /*GMBOLUS*/
         ENDDO
        ENDDO
       ENDDO
       DO J=1,JE
        DO I=1,IE
         FLUM(I,J)=ZERO
         PEM(I,J)=ZERO
         SICTRU(I,J)=ZERO
         SICTRV(I,J)=ZERO
         SUM_EMINPO(I,J)=ZERO
         SUM_ZO(I,J)=ZERO
         SUM_FLUM(I,J)=ZERO
         SUM_PEM(I,J)=ZERO
         SUM_SICTHO(I,J)=ZERO
         SUM_SICOMO(I,J)=ZERO
         SUM_SICUO(I,J)=ZERO
         SUM_SICVE(I,J)=ZERO
         SUM_SICTRU(I,J)=ZERO
         SUM_SICTRV(I,J)=ZERO
         SUM_SICSNO(I,J)=ZERO
         SUM_PRECH(I,J)=ZERO
         SUM_RIVRUN(I,J)=ZERO
         SUM_QSWO(I,J)=ZERO
         SUM_QLWO(I,J)=ZERO
         SUM_QLAO(I,J)=ZERO
         SUM_QSEO(I,J)=ZERO
        ENDDO
       ENDDO
#endif /*MEAN*/

#ifdef NON_HYDROSTATIC
#ifdef PARTIALCELL
      DO J=1,JE
      DO I=1,IE
      DO K=1,KEP
        DZPNH(I,J,K) = 0.0
      ENDDO
      ENDDO
      ENDDO
#endif /*PARTIALCELL*/

      DO J=1,JE
      DO I=I_start,IE
      DO K=1,KE
        UOO_1(I,J,K)=0.0
        GU(I,J,K)=0.0
      ENDDO
      ENDDO
      ENDDO

      DO J=J_start,JE
      DO I=1,IE
      DO K=1,KE
        VOE_1(I,J,K)=0.0
        GV(I,J,K)=0.0
      ENDDO
      ENDDO
      ENDDO

      DO J=1,JE
      DO I=1,IE
      DO K=1,KEP
        WNO(I,J,K)=0.0
        WN1O(I,J,K)=0.0
        WNO_1(I,J,K)=0.0
        GW(I,J,K)=0.0        
      ENDDO
      ENDDO
      ENDDO

      DO J=1,JE
      DO I=1,IE
      DO K=1,KE
        DIVGU(I,J,K)=0.0
        DIVGV(I,J,K)=0.0
        DIVGW(I,J,K)=0.0
        DIVG(I,J,K)=0.0
        PNH(I,J,K)=0.0               
      ENDDO
      ENDDO
      ENDDO
#ifdef CG2D4NONHY      
      DO J=1,JE
      DO I=I_start,IE
      	BOTGU(I,J)=0.0     	
      ENDDO
      ENDDO
           
      DO J=J_start,JE
      DO I=1,IE
      	BOTGV(I,J)=0.0     	
      ENDDO
      ENDDO      
      
      DO J=1,JE
      DO I=1,IE
      	BOTGW(I,J)=0.0
        DIVG_bot_surf(I,J)=0.0
        PNH_bot(I,J)=0.0
        DIVBOTGU(I,J)=0.0
        DIVBOTGV(I,J)=0.0
        DIVBOTGW(I,J)=0.0      	      	
      ENDDO
      ENDDO 
#endif /*CG2D4NONHY*/ 
           
#endif
!

      RETURN
      END
