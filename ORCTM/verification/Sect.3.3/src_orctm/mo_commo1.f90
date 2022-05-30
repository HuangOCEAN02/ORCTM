      MODULE MO_COMMO1

      USE mo_param1
#ifdef FRAGMENT      
      USE mo_kind
#endif
      IMPLICIT NONE

      INTEGER IMITTEL,IAUFR,IAUFW,JAHR,IKLIM,NNNDT,IMMA,JMMA
      INTEGER ISTART,I3DREST,mmccdt,d_date
      REAL DT,ROCON,dt0jd,rlsa,rlsa1,chovop
!
      REAL ZERO,ONE,TWO,TENM4,FOUR,HALF,FOURTH,EIGHTH
      REAL G,AH,ROCP,ROCD,PI,RADIUS,OMEGA,RELTEM,RELSAL,DH,DTI
      REAL TIPHAS,COPHAS,SIPHAS,GAMMA0
      INTEGER LITER,NDUMMY,N1,M1,N2,M2,N3,M3,N4,M4,KB,KBM,KM,MATR,NMAX, &
     &  IUNIT,IFLAG,LYEARS,LDTYEAR, LYEAR1,LYEAR2,LYEAR,LMONTH,LMONTS,  &
     &  LDTMONTH,LMONT1,LMONT2,MONLEN(12),NDTMONTH,LDAYS,LDAY,LDTDAYC,  &
     &  LDT,NYEARS,NMONTS,LY_END,LY_START,LM_START
!
      REAL ALMZER,STABN,STABO,CONO,CONN,GHO,GHN,DTDECT,CLICN,CLICO,     &
     &  AULAPTS,AULAPUV,AUS,AH00,CAH00,AV0,DV0,WT,CSTABEPS,CAVOCON,     &
     &  CDVOCON
!
      REAL OLDSO
      INTEGER NTIPHA
#ifdef NON_HYDROSTATIC
      LOGICAL FIRST_STEP
#endif
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      REAL, POINTER :: DLXP(:,:),DLYP(:,:),DLXU(:,:),DLYV(:,:),         &
     &   DLXV(:,:),DLYU(:,:),DTDXUO(:,:),DTDXUE(:,:),DTDYO(:,:),        &
     &   DTDXPO(:,:),DPYE(:,:),DPYO(:,:),DDUE(:,:,:),DDUO(:,:,:),       &
     &   DDPO(:,:,:),DPIO(:,:,:),DDPSIO(:,:,:),DWI(:),DEUTE(:,:),       &
     &   DEUTO(:,:),DEUTIE(:,:),DEUTIO(:,:),DEPTO(:,:),SIFTFE(:),       &
     &   SIFTFO(:),COFTFE(:),COFTFO(:)
! RJ: The NAMELIST doesn't work if DZW is a pointer, so
!     arrays depending on KE are declared fixed
      REAL DZ(KE+1),DZW(KE),TIESTU(KE+1),TIESTW(KE+1),DI(KE+1)
!
!:: GEOGRAPHIC POSITION:
      REAL, POINTER :: GILA(:,:),GIPH(:,:),ALAT(:,:)
      REAL, POINTER :: GILA_G(:,:),GIPH_G(:,:)
!
!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      REAL, POINTER :: WETO(:,:,:),AMSUE(:,:,:),AMSUO(:,:,:)
!
!:: ONE-DIMENSIONAL FIELDS
      REAL, POINTER :: DOTRT(:),UPTRT(:),WINTUR(:),DBACKV(:),ABACKV(:), &
     &   TAF(:),SAF(:),PREFF(:),ROREF(:),SKAL(:),B(:),X(:)
!
!:: THREE-DIMENSIONAL FIELDS
      REAL, POINTER :: RHOO(:,:,:),UKO(:,:,:),VKE(:,:,:),THO(:,:,:),    &
     &   SAO(:,:,:),WO(:,:,:),PO(:,:,:),VK1E(:,:,:),UK1O(:,:,:),        &
     &   VOE(:,:,:),UOO(:,:,:),T1O(:,:,:),S1O(:,:,:),DVO(:,:,:),        &
     &   AVO(:,:,:),STABIO(:,:,:)
!
!:: TWO-DIMENSIONAL FIELDS
      REAL, POINTER :: ZO(:,:),Z1O(:,:),                                &
     &   V1E(:,:),U1O(:,:),USO(:,:),VSE(:,:),UCOS(:,:),VCOS(:,:),       &
     &   SHELP(:,:),THELP(:,:),RHELP(:,:),TLOW(:,:),SLOW(:,:),          &
     &   TSUP(:,:),SSUP(:,:),B1O(:,:),PXOIN(:,:),PYEIN(:,:),UZO(:,:),   &
     &   VZE(:,:),FTWOU(:,:),FTWOV(:,:)
      INTEGER, POINTER :: NUM(:,:),KBOT(:,:),KCONDEP(:,:),NUM_G(:,:)
!
!:: ICE-MODEL FIELDS
      REAL, POINTER ::                                                  &
     &   SICTHO(:,:),SICOMO(:,:),SICOMU(:,:),SICOMV(:,:),               &
     &   SICOMP(:,:),SICUO(:,:),SICVE(:,:),HIBDELO(:,:),                &
     &   HIBZETO(:,:),HIBETO(:,:),HIBDELE(:,:),HIBZETE(:,:),            &
     &   HIBETE(:,:),TAUWATU(:,:),TAUWATV(:,:),SICSNO(:,:),             &
     &   FWO(:,:),SHO(:,:)
!
!:: FORCING FIELDS
      REAL, POINTER :: TXO(:,:),TYE(:,:),TAFO(:,:),EMINPO(:,:),         &
     &   RELSAO(:,:),RELTHO(:,:),FCLOU(:,:),FSWR(:,:),FU10(:,:),        &
     &   FPREC(:,:),FTDEW(:,:),RIVRUN(:,:),GIRIV(:,:)
!
#ifdef NON_HYDROSTATIC
      REAL,POINTER::WNO(:,:,:),WN1O(:,:,:),UOO_1(:,:,:),VOE_1(:,:,:),   &
     & WNO_1(:,:,:),GU(:,:,:),GV(:,:,:),GW(:,:,:),DIVGU(:,:,:),         &
     & DIVGV(:,:,:),DIVGW(:,:,:),DIVG(:,:,:),PNH(:,:,:),DZPNH(:,:,:)            
          
#ifdef CG2D4NONHY      
      REAL,POINTER::PNH_bot(:,:),DIVG_bot_surf(:,:),BOTGU(:,:),         &
     & BOTGV(:,:),BOTGW(:,:),DIVBOTGU(:,:),DIVBOTGV(:,:),DIVBOTGW(:,:)
#endif      
#endif
!
#ifdef GMBOLUS
      REAL, POINTER :: WGO(:,:,:),BOLX(:,:),BOLY(:,:)
#endif
!
#ifdef FRAGMENT
      REAL(Kind=sp),ALLOCATABLE :: DUMMY_P(:,:),DUMMY_U(:,:),         &
     & DUMMY_V(:,:)
#endif
!
      CONTAINS

      SUBROUTINE alloc_mem_commo1
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      ALLOCATE( DLXP(IE,JE),DLYP(IE,JE),DLXU(I_start:IE,JE),            &
     &  DLYV(IE,J_start:JE),DLXV(IE,J_start:JE),DLYU(I_start:IE,JE),    &
     &  DTDXUO(I_start:IE,JE),DTDXUE(IE,J_start:JE),DTDYO(IE,JE),       &
     &  DTDXPO(IE,JE),DPYE(IE,J_start:JE),DPYO(I_start:IE,JE),          &
     &  DDUE(IE,J_start:JE,KE),DDUO(I_start:IE,JE,KE),DDPO(IE,JE,KE),   &
     &  DPIO(IE,JE,KE),DDPSIO(IE,JE,KE),DWI(KE),DEUTE(IE,J_start:JE),   &
     &  DEUTO(I_start:IE,JE),DEUTIE(IE,J_start:JE),                     &
     &  DEUTIO(I_start:IE,JE),DEPTO(IE,JE),                             &
     &  SIFTFE(JE),SIFTFO(JE),COFTFE(JE),COFTFO(JE))
!
!:: GEOGRAPHIC POSITION:
      ALLOCATE(GILA(ITO,JTO),GIPH(ITO,JTO),ALAT(IE,JE))
      ALLOCATE(GILA_G(2*IE_G+1,2*JE_G+1),GIPH_G(2*IE_G+1,2*JE_G+1))
!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      ALLOCATE(AMSUE(IE,J_start:JE,KE),AMSUO(I_start:IE,JE,KE),         &
     &   WETO(IE,JE,KE))
!
!:: ONE-DIMENSIONAL FIELDS
      ALLOCATE( DOTRT(KE),UPTRT(KE),WINTUR(KEP),DBACKV(KEP),            &
     &   ABACKV(KEP),TAF(KE),SAF(KE),PREFF(KE),ROREF(KE),SKAL(ILL),     &
     &   B(ILT),X(ILT))
!
!:: THREE-DIMENSIONAL FIELDS
      ALLOCATE( RHOO(IE,JE,KE),UKO(I_start:IE,JE,KE),                   &
     &   VKE(IE,J_start:JE,KE),THO(IE,JE,KE),SAO(IE,JE,KE),             &
     &   WO(IE,JE,KEP),PO(IE,JE,KE),VK1E(IE,J_start:JE,KE),             &
     &   UK1O(I_start:IE,JE,KE),VOE(IE,J_start:JE,KE),                  &
     &   UOO(I_start:IE,JE,KE),T1O(IE,JE,KE),S1O(IE,JE,KE),             &
     &   DVO(IE,JE,KEP),AVO(IE,JE,KEP),STABIO(IE,JE,KE))

#ifdef NON_HYDROSTATIC
      ALLOCATE(PNH(IE,JE,KE),WNO(IE,JE,KEP),WN1O(IE,JE,KEP),            &
     &  UOO_1(I_start:IE,JE,KE),GU(I_start:IE,JE,KE),DIVGU(IE,JE,KE),   &
     &  VOE_1(IE,J_start:JE,KE),GV(IE,J_start:JE,KE),DIVGV(IE,JE,KE),   &
     &  WNO_1(IE,JE,KEP),GW(IE,JE,KEP),DIVGW(IE,JE,KE),DIVG(IE,JE,KE),  &
     &  DZPNH(IE,JE,KEP))      
#ifdef CG2D4NONHY       
     ALLOCATE(DIVG_bot_surf(IE,JE),PNH_bot(IE,JE),BOTGU(I_start:IE,JE), &
     &  BOTGV(IE,J_start:JE),BOTGW(IE,JE),DIVBOTGU(IE,JE),              &
     &  DIVBOTGV(IE,JE),DIVBOTGW(IE,JE))
#endif /*CG2D4NONHY*/      
#endif /*NON_HYDROSTATIC*/
!
!:: TWO-DIMENSIONAL FIELDS
      ALLOCATE( ZO(IE,JE),Z1O(IE,JE),                                   &
     & V1E(IE,J_start:JE),U1O(I_start:IE,JE),USO(I_start:IE,JE),        &
     & VSE(IE,J_start:JE),UCOS(I_start:IE,JE),VCOS(IE,J_start:JE),      &
     & SHELP(IE,JE),THELP(IE,JE),RHELP(IE,JE),TLOW(IE,JE),SLOW(IE,JE),  &
     & TSUP(IE,JE),SSUP(IE,JE),B1O(IE,JE),PXOIN(I_start:IE,JE),         &
     & PYEIN(IE,J_start:JE),UZO(I_start:IE,JE),VZE(IE,J_start:JE),      &
     & FTWOU(I_start:IE,JE),FTWOV(IE,J_start:JE) )
      ALLOCATE(NUM(IE,JE),KBOT(IE,JE),KCONDEP(IE,JE),NUM_G(IE_G,JE_G))
!
!:: ICE-MODEL FIELDS
       ALLOCATE(                                                        &
     &   SICTHO(IE,JE),SICOMO(IE,JE),SICOMU(IE,JE),SICOMV(IE,JE),       &
     &   SICOMP(IE,JE),SICUO(IE,JE),SICVE(IE,JE),HIBDELO(IE,JE),        &
     &   HIBZETO(IE,JE),HIBETO(IE,JE),HIBDELE(IE,JE),HIBZETE(IE,JE),    &
     &   HIBETE(IE,JE),TAUWATU(IE,JE),TAUWATV(IE,JE),SICSNO(IE,JE),     &
     &   FWO(IE,JE),SHO(IE,JE))
!
!:: FORCING FIELDS
      ALLOCATE( TXO(I_start:IE,JE),TYE(IE,J_start:JE),TAFO(IE,JE),      &
     &    EMINPO(IE,JE),RELSAO(IE,JE),RELTHO(IE,JE),FCLOU(IE,JE),       &
     &    FSWR(IE,JE),FU10(IE,JE),FPREC(IE,JE),FTDEW(IE,JE),            &
     &    RIVRUN(IE,JE),GIRIV(IE,JE))
!
#ifdef GMBOLUS
      ALLOCATE(WGO(IE,JE,KEP),BOLX(IE,JE),BOLY(IE,JE))
#endif
!
#ifdef FRAGMENT
      ALLOCATE(DUMMY_P(IE,JE),DUMMY_U(I_start:IE,JE),      &
     &    DUMMY_V(IE,J_start:JE))
#endif
!
!<===END COMMO1
      END SUBROUTINE alloc_mem_commo1

      END MODULE MO_COMMO1
