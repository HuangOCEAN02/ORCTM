      MODULE MO_OBCS_NEST
!=======================================================================
!   19.03.2023 update for vertical velocity open boundaryconditions                    !
!   17.04.2023 update for Flather and Orlanski radiation conditions                     !
!   20.05.2023 update for Tidal and subtidal open boundary conditions               !
!========================================== Hao Huang====================
!                                                                      !
!  These subroutines set Open Boundary conditions for tidal and subtidal nesting simulation    !
!
!=======================================================================
!    -----------------------ORCTM-------------------------
!     CALL alloc_mem_obcs_nest
!     #ifdef OBCS_TST_NEST
!       CALL alloc_mem_obcs_tst
!     #endif
!
!     CALL INI_OBCS_NEST
!     #ifdef OBCS_TST_NEST
!        CALL INI_OBCS_TST
!     #endif
!
!    #ifdef OBCS_TST_NEST
!      CALL FORC_TST_TIDES(LDTRUN)
!   #endif /*OBCS_TST_NEST*/
!
!     CALL FORC_OBCS_ETA
!         CALL UPGRADE_FLA_BT_ZO
!        #ifdef OBC_Subtide
!           CALL UPGRADE_FLA_SUBT_ZO
!         #endif
!
!     CALL FORC_OBCS_UV
!         CALL UPGRADE_FLA_BT_UV
!        #ifdef OBC_Subtide
!           CALL UPGRADE_ORL_SUBT_UV
!         #endif
!
!    #ifdef OBCS_TST_NEST
!     CALL SET_TIDES
!    #endif
!
!    ---------------------OCVTRO--------------------------
!     CALL OBCS_CHAPMAN
!   #ifdef OBCS_TST_NEST
!       CALL OBCS_TST_UVO(I_step)
!    #else  ! Flather boundary for barotropic velocity
!       CALL OBCS_FLATHER
!    #endif
!    -----------------------ORCTM-------------------------
! #if defined ORLANSKI && defined OBCS_UV_FLATHER 
! #ifdef OBCS_TST_NEST
!  Tidal and subtidal OBC basd on Liu and Gan 2016
!     CALL OBCS_TST_UVK(UOO,VOE,LDTRUN)
!     CALL OBCS_TST_TS(THO,SAO,LDTRUN)
! #else
!   Nudging or Radiation OBC based on ORLANSKI (1976)
!     CALL ORLANSKI_UV(UOO,VOE,LDTRUN)
!     CALL ORLANSKI_TS(THO,SAO,LDTRUN)
! #endif
! #endif
!=======================================================================
      USE mo_commo1
      USE mo_param1
      USE mo_mpi
      USE mo_parallel
      USE mo_units
      USE mo_obcs
!!
      REAL, POINTER :: OBWW(:,:),OBEW(:,:),OBNW(:,:),OBSW(:,:)
      REAL, POINTER :: obc_zwt(:),obc_zet(:),obc_znt(:),obc_zst(:)
!
      REAL, POINTER :: obc_zws(:),obc_zes(:),obc_zns(:), &
      &            obc_zss(:) ! subtidal sea level
      REAL, POINTER :: obcf_uw(:,:),obcf_ue(:,:),obcf_un(:,:),          &
     &                 obcf_us(:,:)
      REAL, POINTER :: obcf_vw(:,:),obcf_ve(:,:),obcf_vn(:,:),          &
     &                 obcf_vs(:,:)
      REAL, POINTER :: OBWUsubtide(:,:),OBEUsubtide(:,:), &
      &                              OBNUsubtide(:,:),OBSUsubtide(:,:), &
      &                              OBWVsubtide(:,:),OBEVsubtide(:,:), &
      &                              OBNVsubtide(:,:),OBSVsubtide(:,:), &
      &                              OBWZsubtide(:)  ,OBEZsubtide(:),   &
      &                              OBNZsubtide(:)  ,OBSZsubtide(:)

#ifdef OBCS_UV_FLATHER
      REAL, POINTER :: ZOOLD(:,:),USOOLD(:,:),VSEOLD(:,:),zetaT(:,:),   &
     &                USOT(:,:),VSET(:,:),ZONEW(:,:)

      REAL, POINTER :: zeta_subT(:,:)
      REAL, POINTER :: obc_ubarw(:),obc_ubare(:),obc_ubarn(:),obc_ubars(:)
      REAL, POINTER :: obc_vbarw(:),obc_vbare(:),obc_vbarn(:),obc_vbars(:)

#endif /*OBCS_UV_FLATHER*/

#ifdef ORLANSKI
      REAL, POINTER :: gradu_wb(:,:),gradu_eb(:,:),gradu_sb(:,:),gradu_nb(:,:),  &
     &                 gradv_wb(:,:),gradv_eb(:,:),gradv_sb(:,:),gradv_nb(:,:),   &
     &                 grad_wb(:,:),grad_eb(:,:), grad_sb(:,:),grad_nb(:,:)

      REAL dUdt, dVdt, dUdx, dVdx, dUde, dVde, dTdt, dTdx ,dTde, cff,Cx &
     &, Ce, dSdt, dSdx, dSde, elevation_tiaozheng, cff2, tau, phi
#endif  /*ORLANSKI*/

#ifdef OBCS_W_FORC
      INTEGER,PARAMETER::                                                   &
     & IO_IN_OBCWNHW=919,IO_IN_OBCWNHE=929,IO_IN_OBCWNHN=939,IO_IN_OBCWNHS=949
#endif /*OBCS_W_FORC*/

#ifdef OBCS_TST_NEST
     INTEGER,PARAMETER::                                               &
     &        IO_IN_TSTUAW=949,IO_IN_TSTUGW=950,  &
     &         IO_IN_TSTUAE=951,IO_IN_TSTUGE=952,   &
     &         IO_IN_TSTUAN=953,IO_IN_TSTUGN=954,   &
     &         IO_IN_TSTUAS=955,IO_IN_TSTUGS=956,   &
     &         IO_IN_TSTVAW=957,IO_IN_TSTVGW=958,  &
     &         IO_IN_TSTVAE=959,IO_IN_TSTVGE=960,   &
     &         IO_IN_TSTVAN=961,IO_IN_TSTVGN=962,   &
     &         IO_IN_TSTVAS=963,IO_IN_TSTVGS=964,   &
     &         IO_IN_TSTZAW=965,IO_IN_TSTZGW=966,  &
     &         IO_IN_TSTZAE=967,IO_IN_TSTZGE=968,   &
     &         IO_IN_TSTZAN=969,IO_IN_TSTZGN=970,   &
     &         IO_IN_TSTZAS=971,IO_IN_TSTZGS=972

     REAL, POINTER :: OBTSTWUam(:,:,:),OBTSTWUph(:,:,:),  &
     &                      OBTSTWVam(:,:,:),OBTSTWVph(:,:,:),  &
     &                      OBTSTWZam(:,:,:),OBTSTWZph(:,:,:),  &
     &                      OBTSTEUam(:,:,:),OBTSTEUph(:,:,:),  &
     &                      OBTSTEVam(:,:,:),OBTSTEVph(:,:,:),  &
     &                      OBTSTEZam(:,:,:),OBTSTEZph(:,:,:),  &
     &                      OBTSTNUam(:,:,:),OBTSTNUph(:,:,:),  &
     &                      OBTSTNVam(:,:,:),OBTSTNVph(:,:,:),  &
     &                      OBTSTNZam(:,:,:),OBTSTNZph(:,:,:),  &
     &                      OBTSTSUam(:,:,:),OBTSTSUph(:,:,:),  &
     &                      OBTSTSVam(:,:,:),OBTSTSVph(:,:,:),  &
     &                      OBTSTSZam(:,:,:),OBTSTSZph(:,:,:)

     REAL, POINTER :: Utidew(:,:),Etidew(:,:),Vtidew(:,:),  &
     &                      Utidee(:,:),Etidee(:,:),Vtidee(:,:),         &
     &                      Utiden(:,:),Etiden(:,:),Vtiden(:,:),          &
     &                      Utides(:,:),Etides(:,:),Vtides(:,:)

! Barotropic Variables
      REAL, POINTER :: USOR_wb0(:,:),USOR_wb1(:,:),      &
     &                 USOR_eb0(:,:),USOR_eb1(:,:),                     &
     &                 VSER_nb0(:,:),VSER_nb1(:,:),                     &
     &                 VSER_sb0(:,:),VSER_sb1(:,:)
      REAL, POINTER :: obc_ubarNw(:),obc_zsubNw(:),obc_zTNw(:),&
     &                 obc_ubarNe(:),obc_zsubNe(:),obc_zTNe(:),     &
     &                 obc_vbarNs(:),obc_zsubNs(:),obc_zTNs(:),      &
     &                 obc_vbarNn(:),obc_zsubNn(:),obc_zTNn(:)
! Baroclinc Variables
      REAL, POINTER :: uko_wb0(:,:,:),uko_wb1(:,:,:),      &
     &                 uko_eb0(:,:,:),uko_eb1(:,:,:),                     &
     &                 uko_sb0(:,:,:),uko_sb1(:,:,:),                      &
     &                 uko_nb0(:,:,:),uko_nb1(:,:,:),                     &
     &                 vke_wb0(:,:,:),vke_wb1(:,:,:),                    &
     &                 vke_eb0(:,:,:),vke_eb1(:,:,:),                     &
     &                 vke_nb0(:,:,:),vke_nb1(:,:,:),                    &
     &                 vke_sb0(:,:,:),vke_sb1(:,:,:)
      REAL, POINTER  :: obc_uNw(:,:), obc_uTNw(:,:),        &
     &                 obc_vNw(:,:), obc_vTNw(:,:),                    &
     &                 obc_uNe(:,:),obc_uTNe(:,:),                     &
     &                 obc_vNe(:,:), obc_vTNe(:,:),                     &
     &                 obc_uNs(:,:),obc_uTNs(:,:),                      &
     &                 obc_vNs(:,:), obc_vTNs(:,:),                      &
     &                 obc_uNn(:,:),obc_uTNn(:,:),                       &
     &                 obc_vNn(:,:), obc_vTNn(:,:)
      REAL, POINTER  :: cli_tNw(:,:),cli_tNe(:,:),         &
     &                 cli_tNs(:,:),cli_tNn(:,:),                        &
      &                cli_sNw(:,:),cli_sNe(:,:),                      &
     &                 cli_sNs(:,:),cli_sNn(:,:)

#endif  /*OBCS_TST_NEST*/

     CONTAINS

! allocate memory for OBCS module
      SUBROUTINE alloc_mem_obcs_nest

      allocate( OBWW(je_g,kep),OBEW(je_g,kep),OBNW(ie_g,kep),OBSW(ie_g,kep)  )
      allocate( obc_zwt(je),obc_zet(je),obc_znt(ie),obc_zst(ie) )

       allocate(obc_zws(je),obc_zes(je),obc_zns(ie),obc_zss(ie),  &
     &      obcf_uw(je,ke),obcf_ue(je,ke),obcf_un(I_start:ie,ke),     &
     &      obcf_us(I_start:ie,ke), obcf_vw(J_start:je,ke),                &
     &      obcf_ve(J_start:je,ke), obcf_vn(ie,ke), obcf_vs(ie,ke)  )

       allocate( OBWUsubtide(je_g,      ke), OBEUsubtide(je_g,      ke), &
     &          OBNUsubtide(0:ie_g,    ke), OBSUsubtide(0:ie_g,    ke), &
     &          OBWVsubtide(0:je_g,    ke), OBEVsubtide(0:je_g,    ke), &
     &          OBNVsubtide(ie_g,      ke), OBSVsubtide(ie_g,      ke), &
     &          OBWZsubtide(je_g)         , OBEZsubtide(je_g)         , &
     &          OBNZsubtide(ie_g)         , OBSZsubtide(ie_g) )


#ifdef OBCS_UV_FLATHER
      allocate(ZOOLD(IE,JE),USOOLD(I_start:IE,JE),VSEOLD(IE,J_start:JE), &
     &        zetaT(IE,JE),USOT(I_start:IE,JE),VSET(IE,J_start:JE)  )
      allocate(zeta_subT(IE,JE),ZONEW(IE,JE))

      allocate( obc_ubarw(je),obc_ubare(je),       &
     &         obc_ubarn(I_start:ie),obc_ubars(I_start:ie),               &
     &         obc_vbarw(J_start:je),obc_vbare(J_start:je),             &
     &         obc_vbarn(ie),obc_vbars(ie) )

#endif /*OBCS_UV_FLATHER*/

#ifdef ORLANSKI
! The index changes for dU in tangential direction!
      allocate( gradu_wb(0:2,0:je),   gradu_eb(ie1:ie,0:je),               &
     &         gradu_sb(0:ie,je1:je),gradu_nb(0:ie,1:3)      )
! The index changes for dV in tangential direction !
      allocate( gradv_wb(1:3,1:je),   gradv_eb(ie1:ie,1:je),               &
     &         gradv_sb(1:ie+1,je1:je),gradv_nb(1:ie+1,0:2)      )
! The index changes for dP in tangential direction!
      allocate( grad_wb(0:2,1:je),   grad_eb(ie1:ie,1:je),               &
     &         grad_sb(0:ie,je1:je),grad_nb(0:ie,0:2)      )
#endif /*ORLANSKI*/

      OBWW=0.0;OBEW=0.0;OBNW=0.0;OBSW=0.0
      obc_zwt(:)=0.0; obc_zet(:)=0.0;obc_znt(:)=0.0;obc_zst(:)=0.0
      obc_zws(:)=0.0;obc_zes(:)=0.0;obc_zns(:)=0.0;obc_zss(:)=0.0
      obcf_uw=0.0;obcf_ue=0.0;obcf_un=0.0;obcf_us=0.0
      obcf_vw=0.0;obcf_ve=0.0;obcf_vn=0.0;obcf_vs=0.0
      OBWUsubtide(:,:)=0.0; OBEUsubtide(:,:)=0.0
      OBNUsubtide(:,:)=0.0; OBSUsubtide(:,:)=0.0
      OBWVsubtide(:,:)=0.0; OBEVsubtide(:,:)=0.0
      OBNVsubtide(:,:)=0.0; OBSVsubtide(:,:)=0.0
      OBWZsubtide(:)=0.0;   OBEZsubtide(:)=0.0
      OBNZsubtide(:)=0.0;   OBSZsubtide(:)=0.0

#ifdef OBCS_UV_FLATHER
      ZOOLD=0.0; zetaT=0.0 ; USOOLD=0.0; VSEOLD=0.0
      USOT=0.0; VSET=0.0
      zeta_subT=0.0;ZONEW=0.0
      obc_ubarw=0.0;obc_ubare=0.0;obc_ubarn=0.0;obc_ubars=0.0
      obc_vbarw=0.0;obc_vbare=0.0;obc_vbarn=0.0;obc_vbars=0.0
#endif /*OBCS_UV_FLATHER*/

#ifdef ORLANSKI
      gradu_wb=0.0;gradu_eb=0.0;gradu_sb=0.0;gradu_nb=0.0
      gradv_wb=0.0;gradv_eb=0.0;gradv_sb=0.0;gradv_nb=0.0
      grad_wb=0.0;grad_eb=0.0;grad_sb=0.0;grad_nb=0.0
      dUdt=0.0; dVdt=0.0; dUdx=0.0; dVdx=0.0; dUde=0.0; dVde=0.0; dTdt=0.0
      dTdx=0.0; dTde=0.0; cff=0.0;; Cx=0.0 ;  Ce=0.0; dSdt=0.0; dSdx=0.0
      dSde=0.0; elevation_tiaozheng=0.0; cff2=0.0; tau=0.0; phi=0.0
#endif  /*ORLANSKI*/

      END  SUBROUTINE alloc_mem_obcs_nest

     SUBROUTINE alloc_mem_obcs_tst
#ifdef OBCS_TST_NEST
     allocate( OBTSTWZam(1:je_g,1:tidalnumber,3),OBTSTEZam(1:je_g,1:tidalnumber,3), &
     &         OBTSTNZam(1:ie_g,1:tidalnumber,3),OBTSTSZam(1:ie_g,1:tidalnumber,3), &
     &         OBTSTWZph(1:je_g,1:tidalnumber,3),OBTSTEZph(1:je_g,1:tidalnumber,3), &
     &         OBTSTNZph(1:ie_g,1:tidalnumber,3),OBTSTSZph(1:ie_g,1:tidalnumber,3), &

     &         OBTSTWUam(1:je_g,1:tidalnumber,3),OBTSTEUam(1:je_g,1:tidalnumber,3), &
     &         OBTSTNUam(0:ie_g,1:tidalnumber,3),OBTSTSUam(0:ie_g,1:tidalnumber,3), &
     &         OBTSTWUph(1:je_g,1:tidalnumber,3),OBTSTEUph(1:je_g,1:tidalnumber,3), &
     &         OBTSTNUph(0:ie_g,1:tidalnumber,3),OBTSTSUph(0:ie_g,1:tidalnumber,3), &

     &         OBTSTWVam(0:je_g,1:tidalnumber,3),OBTSTEVam(0:je_g,1:tidalnumber,3), &
     &         OBTSTNVam(1:ie_g,1:tidalnumber,3),OBTSTSVam(1:ie_g,1:tidalnumber,3), &
     &         OBTSTWVph(0:je_g,1:tidalnumber,3),OBTSTEVph(0:je_g,1:tidalnumber,3), &
     &         OBTSTNVph(1:ie_g,1:tidalnumber,3),OBTSTSVph(1:ie_g,1:tidalnumber,3))

      ALLOCATE(Utidew(1:3,1:je),Vtidew(1:3,J_start:je),Etidew(1:3,1:je),  &
     & Utidee(ie-3:ie1,1:je),Vtidee(ie2:ie,J_start:je),Etidee(ie2:ie,1:je),    &
     & Utiden(I_start:ie,1:3),Vtiden(1:ie,1:3),Etiden(1:ie,1:3),                   &
     & Utides(I_start:ie,je2:je),Vtides(1:ie,je-3:je1),Etides(1:ie,je2:je))

! Barotropic Variables
      ALLOCATE(USOR_wb0(1:3,je),USOR_wb1(1:3,je),      &
     &                 USOR_eb0(ie-3:ie1,je),USOR_eb1(ie-3:ie1,je),          &
     &                 VSER_nb0(ie,1:3),VSER_nb1(ie,1:3),                     &
     &                 VSER_sb0(ie,je-3:je1),VSER_sb1(ie,je-3:je1))
     ALLOCATE(obc_ubarNw(je),obc_zsubNw(je),obc_zTNw(je),&
     &                 obc_ubarNe(je),obc_zsubNe(je),obc_zTNe(je),     &
     &                 obc_vbarNs(ie),obc_zsubNs(ie),obc_zTNs(ie),      &
     &                 obc_vbarNn(ie),obc_zsubNn(ie),obc_zTNn(ie))
! Baroclinic Variables
      ALLOCATE(uko_wb0(1:3,je,ke),uko_wb1(1:3,je,ke),      &
     &                 uko_eb0(ie-3:ie1,je,ke),uko_eb1(ie-3:ie1,je,ke),   &
     &                 uko_sb0(i_start:ie,je2:je,ke),uko_sb1(i_start:ie,je2:je,ke), &
     &                 uko_nb0(i_start:ie,1:3,ke),uko_nb1(i_start:ie,1:3,ke),       &
     &                 vke_wb0(1:3,J_start:je,ke),vke_wb1(1:3,J_start:je,ke),      &
     &                 vke_eb0(ie2:ie,J_start:je,ke),vke_eb1(ie2:ie,J_start:je,ke),&
     &                 vke_nb0(ie,1:3,ke),vke_nb1(ie,1:3,ke),                    &
     &                 vke_sb0(ie,je-3:je1,ke),vke_sb1(ie,je-3:je1,ke))
      ALLOCATE(obc_uNw(je,ke), obc_uTNw(1:3,je),        &
     &                 obc_vNw(J_start:je,ke), obc_vTNw(1:3,J_start:je),    &
     &                 obc_uNe(je,ke),obc_uTNe(ie-3:ie1,je),                       &
     &                 obc_vNe(J_start:je,ke), obc_vTNe(ie2:ie,J_start:je),  &
     &                 obc_uNs(I_start:ie,ke),obc_uTNs(I_start:ie,je2:je),     &
     &                 obc_vNs(ie,ke), obc_vTNs(1:ie,je-3:je-1),          &
     &                 obc_uNn(I_start:ie,ke),obc_uTNn(I_start:ie,1:3),       &
     &                 obc_vNn(ie,ke), obc_vTNn(1:ie,1:3))   ! obc_vTNs obc_vTNn bug has found
      ALLOCATE(cli_tNw(je,ke),cli_tNe(je,ke),         &
     &                 cli_tNs(ie,ke),cli_tNn(ie,ke),           &
     &                 cli_sNw(je,ke),cli_sNe(je,ke),        &
     &                 cli_sNs(ie,ke),cli_sNn(ie,ke))

      OBTSTWZam=0.0;OBTSTEZam=0.0;OBTSTNZam=0.0;OBTSTSZam=0.0
      OBTSTWZph=0.0;OBTSTEZph=0.0;OBTSTNZph=0.0;OBTSTSZph=0.0
      OBTSTWUam=0.0;OBTSTEUam=0.0;OBTSTNUam=0.0;OBTSTSUam=0.0
      OBTSTWUph=0.0;OBTSTEUph=0.0;OBTSTNUph=0.0;OBTSTSUph=0.0
      OBTSTWVam=0.0;OBTSTEVam=0.0;OBTSTNVam=0.0;OBTSTSVam=0.0
      OBTSTWVph=0.0;OBTSTEVph=0.0;OBTSTNVph=0.0;OBTSTSVph=0.0

      Utidew=0.0;Vtidew=0.0;Etidew=0.0;Utidee=0.0;Vtidee=0.0;Etidee=0.0
      Utiden=0.0;Vtiden=0.0;Etiden=0.0;Utides=0.0;Vtides=0.0;Etides=0.0
! Barotropic Variables
      USOR_wb0=0.0;USOR_wb1=0.0;USOR_eb0=0.0;USOR_eb1=0.0
      VSER_nb0=0.0;VSER_nb1=0.0;VSER_sb0=0.0;VSER_sb1=0.0

      obc_ubarNw=0.0;obc_zsubNw=0.0;obc_zTNw=0.0
      obc_ubarNe=0.0;obc_zsubNe=0.0;obc_zTNe=0.0
      obc_vbarNn=0.0;obc_zsubNn=0.0;obc_zTNn=0.0
      obc_vbarNs=0.0;obc_zsubNs=0.0;obc_zTNs=0.0
! Baroclinic Variables
      uko_wb0=0.0;uko_wb1=0.0;vke_wb0=0.0;vke_wb1=0.0
      uko_eb0=0.0;uko_eb1=0.0;vke_eb0=0.0;vke_eb1=0.0
      uko_sb0=0.0;uko_sb1=0.0;vke_sb0=0.0;vke_sb1=0.0
      uko_nb0=0.0;uko_nb1=0.0;vke_nb0=0.0;vke_nb1=0.0

      obc_uNw=0.0;obc_uNe=0.0;obc_uNn=0.0;obc_uNs=0.0
      obc_vNw=0.0;obc_vNe=0.0;obc_vNn=0.0;obc_vNs=0.0
      obc_uTNw=0.0;obc_uTNe=0.0;obc_uTNn=0.0;obc_uTNs=0.0
      obc_vTNw=0.0;obc_vTNe=0.0;obc_vTNn=0.0;obc_vTNs=0.0

      cli_tNw=0.0;cli_sNw=0.0;cli_tNe=0.0;cli_sNe=0.0
      cli_tNn=0.0;cli_sNn=0.0;cli_tNs=0.0;cli_sNs=0.0
#endif  /*OBCS_TST_NEST*/
      END SUBROUTINE alloc_mem_obcs_tst

      SUBROUTINE INI_OBCS_NEST

#if defined NON_HYDROSTATIC && defined OBCS_W_FORC
      IF(p_pe==p_io) THEN

       IF( West ) THEN
        OPEN(IO_IN_OBCWNHW,FILE='./OBCW/OB_West_W.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF

       IF( East ) THEN
        OPEN(IO_IN_OBCWNHE,FILE='./OBCW/OB_East_W.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF

       IF( North ) THEN
        OPEN(IO_IN_OBCWNHN,FILE='./OBCW/OB_North_W.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF

       IF( South ) THEN
        OPEN(IO_IN_OBCWNHS,FILE='./OBCW/OB_South_W.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF

       ENDIF
#endif

      END SUBROUTINE INI_OBCS_NEST

     SUBROUTINE INI_OBCS_TST
#ifdef OBCS_TST_NEST
      IF(p_pe==p_io) THEN
       IF( West ) THEN
        OPEN(IO_IN_TSTUAW,FILE='./OB_tst/OBTST_West_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTUGW,FILE='./OB_tst/OBTST_West_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTVAW,FILE='./OB_tst/OBTST_West_v_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTVGW,FILE='./OB_tst/OBTST_West_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTZAW,FILE='./OB_tst/OBTST_West_z_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTZGW,FILE='./OB_tst/OBTST_West_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        READ(IO_IN_TSTUAW)   ((OBTSTWUam(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGW)  ((OBTSTWUph(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAW)   ((OBTSTWUam(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGW)  ((OBTSTWUph(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAW)   ((OBTSTWUam(i,j,3),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGW)  ((OBTSTWUph(i,j,3),i=1,je_g),j=1,tidalnumber)

        READ(IO_IN_TSTVAW)   ((OBTSTWVam(i,j,1),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGW)  ((OBTSTWVph(i,j,1),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAW)   ((OBTSTWVam(i,j,2),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGW)  ((OBTSTWVph(i,j,2),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAW)   ((OBTSTWVam(i,j,3),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGW)  ((OBTSTWVph(i,j,3),i=0,je_g),j=1,tidalnumber)

        READ(IO_IN_TSTZAW)   ((OBTSTWZam(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGW)  ((OBTSTWZph(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAW)   ((OBTSTWZam(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGW)  ((OBTSTWZph(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAW)   ((OBTSTWZam(i,j,3),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGW)  ((OBTSTWZph(i,j,3),i=1,je_g),j=1,tidalnumber)

       ENDIF

       IF( East ) THEN
       OPEN(IO_IN_TSTUAE,FILE='./OB_tst/OBTST_East_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTUGE,FILE='./OB_tst/OBTST_East_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTVAE,FILE='./OB_tst/OBTST_East_v_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTVGE,FILE='./OB_tst/OBTST_East_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTZAE,FILE='./OB_tst/OBTST_East_z_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTZGE,FILE='./OB_tst/OBTST_East_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')


        READ(IO_IN_TSTUAE)   ((OBTSTEUam(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGE)  ((OBTSTEUph(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAE)   ((OBTSTEUam(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGE)  ((OBTSTEUph(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAE)   ((OBTSTEUam(i,j,3),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGE)  ((OBTSTEUph(i,j,3),i=1,je_g),j=1,tidalnumber)

        READ(IO_IN_TSTVAE)   ((OBTSTEVam(i,j,1),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGE)  ((OBTSTEVph(i,j,1),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAE)   ((OBTSTEVam(i,j,2),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGE)  ((OBTSTEVph(i,j,2),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAE)   ((OBTSTEVam(i,j,3),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGE)  ((OBTSTEVph(i,j,3),i=0,je_g),j=1,tidalnumber)

        READ(IO_IN_TSTZAE)   ((OBTSTEZam(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGE)  ((OBTSTEZph(i,j,1),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAE)   ((OBTSTEZam(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGE)  ((OBTSTEZph(i,j,2),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAE)   ((OBTSTEZam(i,j,3),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGE)  ((OBTSTEZph(i,j,3),i=1,je_g),j=1,tidalnumber)

        ENDIF

       IF( South ) THEN
       OPEN(IO_IN_TSTUAS,FILE='./OB_tst/OBTST_South_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTUGS,FILE='./OB_tst/OBTST_South_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTVAS,FILE='./OB_tst/OBTST_South_v_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTVGS,FILE='./OB_tst/OBTST_South_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTZAS,FILE='./OB_tst/OBTST_South_z_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTZGS,FILE='./OB_tst/OBTST_South_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')


        READ(IO_IN_TSTUAS)   ((OBTSTSUam(i,j,1),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGS)  ((OBTSTSUph(i,j,1),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAS)   ((OBTSTSUam(i,j,2),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGS)  ((OBTSTSUph(i,j,2),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAS)   ((OBTSTSUam(i,j,3),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGS)  ((OBTSTSUph(i,j,3),i=0,ie_g),j=1,tidalnumber)

        READ(IO_IN_TSTVAS)   ((OBTSTSVam(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGS)  ((OBTSTSVph(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAS)   ((OBTSTSVam(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGS)  ((OBTSTSVph(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAS)   ((OBTSTSVam(i,j,3),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGS)  ((OBTSTSVph(i,j,3),i=1,ie_g),j=1,tidalnumber)

        READ(IO_IN_TSTZAS)   ((OBTSTSZam(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGS)  ((OBTSTSZph(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAS)   ((OBTSTSZam(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGS)  ((OBTSTSZph(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAS)   ((OBTSTSZam(i,j,3),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGS)  ((OBTSTSZph(i,j,3),i=1,ie_g),j=1,tidalnumber)

       ENDIF

       IF( North ) THEN
        OPEN(IO_IN_TSTUAN,FILE='./OB_tst/OBTST_North_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTUGN,FILE='./OB_tst/OBTST_North_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTVAN,FILE='./OB_tst/OBTST_North_v_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTVGN,FILE='./OB_tst/OBTST_North_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')

        OPEN(IO_IN_TSTZAN,FILE='./OB_tst/OBTST_North_z_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_TSTZGN,FILE='./OB_tst/OBTST_North_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')


        READ(IO_IN_TSTUAN)   ((OBTSTNUam(i,j,1),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGN)  ((OBTSTNUph(i,j,1),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAN)   ((OBTSTNUam(i,j,2),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGN)  ((OBTSTNUph(i,j,2),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUAN)   ((OBTSTNUam(i,j,3),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTUGN)  ((OBTSTNUph(i,j,3),i=0,ie_g),j=1,tidalnumber)

        READ(IO_IN_TSTVAN)   ((OBTSTNVam(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGN)  ((OBTSTNVph(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAN)   ((OBTSTNVam(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGN)  ((OBTSTNVph(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVAN)   ((OBTSTNVam(i,j,3),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTVGN)  ((OBTSTNVph(i,j,3),i=1,ie_g),j=1,tidalnumber)

        READ(IO_IN_TSTZAN)   ((OBTSTNZam(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGN)  ((OBTSTNZph(i,j,1),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAN)   ((OBTSTNZam(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGN)  ((OBTSTNZph(i,j,2),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZAN)   ((OBTSTNZam(i,j,3),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_TSTZGN)  ((OBTSTNZph(i,j,3),i=1,ie_g),j=1,tidalnumber)
       ENDIF
      ENDIF

      IF( West) THEN
       CALL p_bcast(OBTSTWUam,p_io)
       CALL p_bcast(OBTSTWUph,p_io)
       CALL p_bcast(OBTSTWVam,p_io)
       CALL p_bcast(OBTSTWVph,p_io)
       CALL p_bcast(OBTSTWZam,p_io)
       CALL p_bcast(OBTSTWZph,p_io)
      ENDIF

      IF( East) THEN
       CALL p_bcast(OBTSTEUam,p_io)
       CALL p_bcast(OBTSTEUph,p_io)
       CALL p_bcast(OBTSTEVam,p_io)
       CALL p_bcast(OBTSTEVph,p_io)
       CALL p_bcast(OBTSTEZam,p_io)
       CALL p_bcast(OBTSTEZph,p_io)
      ENDIF

      IF( North) THEN
       CALL p_bcast(OBTSTNUam,p_io)
       CALL p_bcast(OBTSTNUph,p_io)
       CALL p_bcast(OBTSTNVam,p_io)
       CALL p_bcast(OBTSTNVph,p_io)
       CALL p_bcast(OBTSTNZam,p_io)
       CALL p_bcast(OBTSTNZph,p_io)
      ENDIF

      IF( South) THEN
       CALL p_bcast(OBTSTSUam,p_io)
       CALL p_bcast(OBTSTSUph,p_io)
       CALL p_bcast(OBTSTSVam,p_io)
       CALL p_bcast(OBTSTSVph,p_io)
       CALL p_bcast(OBTSTSZam,p_io)
       CALL p_bcast(OBTSTSZph,p_io)
      ENDIF
#endif  /*OBCS_TST_NEST*/
     END SUBROUTINE INI_OBCS_TST

! OBC eta forcing
      SUBROUTINE FORC_OBCS_ETA(i_step)

      REAL :: mytime
      INTEGER :: td

      mytime = i_step * DT

#ifdef OBC_TIDE
      ! west OB
      if ( West ) then
       if (have_g_is) then
        obc_zwt(:) = 0.
        DO td=1,tidalnumber
           obc_zwt(:) = obc_zwt(:) + OBWZam(p_joff+1:p_joff+je,td)*weto(1,:,1)* &
             &      COS(2.0 * PI * (myTime-OBWZph(p_joff+1:p_joff+je,td))/tidalPeriod(td))
        ENDDO
       endif
      endif

      ! east OB
      IF ( East ) THEN
       if (have_g_ie) then
        obc_zet(:) = 0.
        DO td=1,tidalnumber
            obc_zet(:) = obc_zet(:) + OBEZam(p_joff+1:p_joff+je,td)*weto(ie,:,1)* &
             &      COS(2.0 * PI * (myTime-OBEZph(p_joff+1:p_joff+je,td))/tidalPeriod(td))
        ENDDO
       endif
      endif


      ! north OB
      IF ( North ) THEN
       if (have_g_js) then
        obc_znt(:) = 0.
        DO td=1,tidalnumber
           obc_znt(:) = obc_znt(:) + OBNZam(p_ioff+1:p_ioff+ie,td)*weto(:,1,1)*  &
             &      COS(2.0 * PI * (myTime-OBNZph(p_ioff+1:p_ioff+ie,td))/tidalPeriod(td))
        ENDDO
       endif
      endif

      ! south OB
      IF ( South ) THEN
       if (have_g_je) then
        obc_zst(:) = 0.
        DO td=1,tidalnumber
           obc_zst(:) = obc_zst(:) + OBSZam(p_ioff+1:p_ioff+ie,td)*weto(:,je,1)*  &
             &      COS(2.0 * PI * (myTime-OBSZph(p_ioff+1:p_ioff+ie,td))/tidalPeriod(td))
        ENDDO
       endif
      endif

#endif /*OBC_TIDE*/

#ifdef OBC_SUBTIDE
      IF ( MOD( i_step,NINT(CLI_STEP/DT) ) .EQ.  1 ) THEN
        IF(p_pe==p_io) THEN
        	  if ( West ) then
               READ(IO_IN_OBCZRW) OBWZsubtide
            endif

            if ( East ) then
               READ(IO_IN_OBCZRE) OBEZsubtide
            endif

            if ( North ) then
               READ(IO_IN_OBCZRN) OBNZsubtide
            endif

            if ( South ) then
               READ(IO_IN_OBCZRS) OBSZsubtide
            endif
       ENDIF

       if ( West ) CALL p_bcast(OBWZsubtide,p_io)
       if ( East ) CALL p_bcast(OBEZsubtide,p_io)
       if ( North) CALL p_bcast(OBNZsubtide,p_io)
       if ( South) CALL p_bcast(OBSZsubtide,p_io)

      ENDIF

      if ( West ) then
        if (have_g_is) then
        	  obc_zws(:)=0
            obc_zws(1:je)=OBWZsubtide(p_joff+1:p_joff+je)*weto(1,:,1)
        endif
      endif

      if ( East ) then
        if (have_g_ie) then
        	  obc_zes(:)=0
            obc_zes(:)=OBEZsubtide(p_joff+1:p_joff+je)*weto(ie,:,1)
        endif
      endif


      if ( North ) then
        if (have_g_js) then
        	  obc_zns(:)=0
            obc_zns(:)=OBNZsubtide(p_ioff+1:p_ioff+ie)*weto(:,1,1)
        endif
      endif


      if ( South ) then
        if (have_g_je) then
            obc_zss(:)=0
            obc_zss(:)=OBSZsubtide(p_ioff+1:p_ioff+ie)*weto(:,je,1)
        endif
      endif
#endif /*OBC_SUBTIDE*/
      END SUBROUTINE FORC_OBCS_ETA


      SUBROUTINE FORC_OBCS_UV(i_step)

      REAL mytime
      INTEGER td

      mytime = i_step * DT

#ifdef OBC_TIDE
      ! west OB
      if ( West .AND. have_g_is) then
         obc_uzw(:)=0.
         obc_vzw(:)=0.
        do td=1,tidalnumber

         obc_uzw(1:je)=obc_uzw(1:je)+ &
          &      OBWUam(p_joff+1:p_joff+je,td)*amsuo(0,1:je,1)* &
          &      COS(2.0 * PI * (myTime-OBWUph(p_joff+1:p_joff+je,td))/tidalPeriod(td))

         obc_vzw(J_start:je)=obc_vzw(J_start:je)+ &
          &      OBWVam(p_joff+J_start:p_joff+je,td)*amsue(1,J_start:je,1)* &
          &      COS(2.0 * PI * (myTime-OBWVph(p_joff+J_start:p_joff+je,td))/tidalPeriod(td))
        enddo
         ! WRITE(IO_STDOUT,*) 'West FORC_OBCS_UV OBWUam(1:je.1:8) : ', OBWUam(p_joff+1:p_joff+5,1:8)
         ! WRITE(IO_STDOUT,*) 'West FORC_OBCS_UV OBWUph(1:je.1:8) : ', OBWUph(p_joff+1:p_joff+5,1:8)
         ! WRITE(IO_STDOUT,*) 'West FORC_OBCS_UV obc_uzw(1:je) : ',  obc_uzw(1:5)
         ! WRITE(IO_STDOUT,*) 'West FORC_OBCS_UV obc_vzw(1:je) : ',  obc_vzw(1:5) 
 !        if (add_tide) then
 !        obc_uzw(:) = obc_uzw(:)*DEUTO(0,1:je)
 !        obc_vzw(:) = obc_vzw(:)*DEUTE(1,J_start:je)
 !      endif
      endif


      ! east OB
      if ( East .AND. have_g_ie) then
         obc_uze(:)=0.
         obc_vze(:)=0.
        do td=1,tidalnumber
          obc_uze(1:je)=obc_uze(1:je)+ &
           &      OBEUam(p_joff+1:p_joff+je,td)*amsuo(ie,1:je,1)* &
           &      COS(2.0 * PI * (myTime-OBEUph(p_joff+1:p_joff+je,td))/tidalPeriod(td))

          obc_vze(J_start:je)=obc_vze(J_start:je)+ &
           &      OBEVam(p_joff+J_start:p_joff+je,td)*amsue(ie,J_start:je,1)* &
           &      COS(2.0 * PI * (myTime-OBEVph(p_joff+J_start:p_joff+je,td))/tidalPeriod(td))
        enddo
 !        if (add_tide) then
 !        obc_uze(:) = obc_uze(:)*DEUTO(ie,1:je)
 !        obc_vze(:) = obc_vze(:)*DEUTE(ie,J_start:je)
 !       endif
      endif


      ! north OB
      if ( North .AND. have_g_js) then
          obc_uzn(:)=0.
          obc_vzn(:)=0.
        do td=1,tidalnumber
          obc_uzn(I_start:ie)=obc_uzn(I_start:ie)+ &
           &      OBNUam(p_ioff+I_start:p_ioff+ie,td)*amsuo(I_start:ie,1,1)* &
           &      COS(2.0 * PI * (myTime-OBNUph(p_ioff+I_start:p_ioff+ie,td))/tidalPeriod(td))

          obc_vzn(1:ie)=obc_vzn(1:ie)+ &
           &      OBNVam(p_ioff+1:p_ioff+ie,td)*amsue(1:ie,0,1)* &
           &      COS(2.0 * PI * (myTime-OBNVph(p_ioff+1:p_ioff+ie,td))/tidalPeriod(td))
        enddo
        
         ! WRITE(IO_STDOUT,*) 'North FORC_OBCS_UV obc_uzn(I_start:ie) : ',  obc_uzn(I_start:5)
         ! WRITE(IO_STDOUT,*) 'North FORC_OBCS_UV obc_vzn(1:ie) : ',  obc_vzn(1:5) 
!        if (add_tide) then
!         obc_uzn(:) = obc_uzn(:)*DEUTO(I_start:ie,1)
!         obc_vzn(:) = obc_vzn(:)*DEUTE(1:ie,0)
!        endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then
         obc_uzs(:)=0.
         obc_vzs(:)=0.
        do td=1,tidalnumber
          obc_uzs(I_start:ie)=obc_uzs(I_start:ie)+ &
           &      OBSUam(p_ioff+I_start:p_ioff+ie,td)*amsuo(I_start:ie,je,1)* &
           &      COS(2.0 * PI * (myTime-OBSUph(p_ioff+I_start:p_ioff+ie,td))/tidalPeriod(td))

          obc_vzs(1:ie)=obc_vzs(1:ie)+ &
           &      OBSVam(p_ioff+1:p_ioff+ie,td)*amsue(1:ie,je,1)* &
           &      COS(2.0 * PI * (myTime-OBSVph(p_ioff+1:p_ioff+ie,td))/tidalPeriod(td))
        enddo
!     if (add_tide) then
!         obc_uzs(:) = obc_uzs(:)*DEUTO(I_start:ie,je)
!         obc_vzs(:) = obc_vzs(:)*DEUTE(1:ie,je)
!      endif
      endif
#endif /*OBC_TIDE*/

#ifdef OBC_SUBTIDE
      IF ( MOD( i_step,NINT(CLI_STEP/DT) ) .EQ.  1 ) THEN
        IF(p_pe==p_io) THEN

          DO K=1,KE
            if ( West ) then
               READ(IO_IN_OBCURW) OBWUsubtide(:,K)
               READ(IO_IN_OBCVRW) OBWVsubtide(:,K)
            endif

            if ( East ) then
               READ(IO_IN_OBCURE) OBEUsubtide(:,K)
               READ(IO_IN_OBCVRE) OBEVsubtide(:,K)
            endif

            if ( North ) then
               READ(IO_IN_OBCURN) OBNUsubtide(:,K)
               READ(IO_IN_OBCVRN) OBNVsubtide(:,K)
            endif

            if ( South ) then
               READ(IO_IN_OBCURS) OBSUsubtide(:,K)
               READ(IO_IN_OBCVRS) OBSVsubtide(:,K)
            endif
          ENDDO
        ENDIF

          if ( West ) then
             CALL p_bcast(OBWUsubtide,p_io)
             CALL p_bcast(OBWVsubtide,p_io)
          endif

          if ( East ) then
             CALL p_bcast(OBEUsubtide,p_io)
             CALL p_bcast(OBEVsubtide,p_io)
          endif

          if ( North ) then
             CALL p_bcast(OBNUsubtide,p_io)
             CALL p_bcast(OBNVsubtide,p_io)
          endif

          if ( South ) then
             CALL p_bcast(OBSUsubtide,p_io)
             CALL p_bcast(OBSVsubtide,p_io)
          endif

       ENDIF

        if ( West .AND. have_g_is) then
          do k=1,ke
            do j=1,je
              obc_usw(j,k)=OBWUsubtide(p_joff+j,k)*amsuo(0,j,k)
            enddo

            do j=J_start,je
              obc_vsw(j,k)=OBWVsubtide(p_joff+j,k)*amsue(1,j,k)
            enddo
          enddo
        endif



        if ( East .AND. have_g_ie) then
          do k=1,ke
            do j=1,je
              obc_use(j,k)=OBEUsubtide(p_joff+j,k)*amsuo(ie,j,k)
            enddo

            do j=J_start,je
              obc_vse(j,k)=OBEVsubtide(p_joff+j,k)*amsue(ie,j,k)
            enddo
          enddo
        endif

        if ( North .AND. have_g_js) then
          do k=1,ke
            do i=I_start,ie
              obc_usn(i,k)=OBNUsubtide(p_ioff+i,k)*amsuo(i,1,k)
            enddo

            do i=1,ie
              obc_vsn(i,k)=OBNVsubtide(p_ioff+i,k)*amsue(i,0,k)
            enddo
          enddo
        endif

        if ( South .AND. have_g_je) then
          do k=1,ke
             do i=I_start,ie
               obc_uss(i,k)=OBSUsubtide(p_ioff+i,k)*amsuo(i,je,k)
             enddo

             do i=1,ie
               obc_vss(i,k)=OBSVsubtide(p_ioff+i,k)*amsue(i,je,k)
             enddo
           enddo
        endif
#endif /*OBC_SUBTIDE*/
      END SUBROUTINE FORC_OBCS_UV

      SUBROUTINE FORC_TST_TIDES(i_step)

      REAL mytime
      INTEGER td
      mytime = i_step * DT
      
#ifdef OBCS_TST_NEST
      ! west OB
      if ( West .AND. have_g_is) then
         Utidew(:,:)=0.
         Vtidew(:,:)=0.
         Etidew(:,:)=0.
        do i = 1,3
        do td=1,tidalnumber
         Utidew(i,1:je)=Utidew(i,1:je)+                                   &
     &      OBTSTWUam(p_joff+1:p_joff+je,td,i)*amsuo(i,1:je,1)*             &
     &      COS(2.0 * PI * (myTime-OBTSTWUph(p_joff+1:p_joff+je,td,i))/tidalPeriod(td))

         Vtidew(i,J_start:je)=Vtidew(i,J_start:je)+                       &
     &      OBTSTWVam(p_joff+J_start:p_joff+je,td,i)*amsue(i,J_start:je,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTWVph(p_joff+J_start:p_joff+je,td,i))/tidalPeriod(td))

         Etidew(i,1:je)=Etidew(i,1:je)+                       &
     &      OBTSTWZam(p_joff+1:p_joff+je,td,i)*weto(i,1:je,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTWZph(p_joff+1:p_joff+je,td,i))/tidalPeriod(td))
        enddo
        enddo
         ! WRITE(IO_STDOUT,*) 'West TST boundary Utidew(1,:) : ',  Utidew(1,1:5)
         ! WRITE(IO_STDOUT,*) 'West TST boundary Vtidew(1,:) : ',  Vtidew(1,1:5)        
         ! WRITE(IO_STDOUT,*) 'West TST boundary Etidew(1,:) : ',  Etidew(1,1:5)        
      endif

      ! east OB
      if ( East .AND. have_g_ie) then
         Utidee(:,:)=0.
         Vtidee(:,:)=0.
         Etidee(:,:)=0.
        do i = ie-3,ie1
        do td=1,tidalnumber
          Utidee(i,1:je)=Utidee(i,1:je)+                                  &
     &      OBTSTEUam(p_joff+1:p_joff+je,td,i-ie+4)*amsuo(i,1:je,1)*           &
     &      COS(2.0 * PI *(myTime-OBTSTEUph(p_joff+1:p_joff+je,td,i-ie+4))/tidalPeriod(td))

          Vtidee(i+1,J_start:je)=Vtidee(i+1,J_start:je)+                      &
     &      OBTSTEVam(p_joff+J_start:p_joff+je,td,i-ie+4)*amsue(i+1,J_start:je,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTEVph(p_joff+J_start:p_joff+je,td,i-ie+4))/tidalPeriod(td))

         Etidee(i+1,1:je)=Etidee(i+1,1:je)+                      &
     &      OBTSTEZam(p_joff+1:p_joff+je,td,i-ie+4)*weto(i+1,1:je,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTEZph(p_joff+1:p_joff+je,td,i-ie+4))/tidalPeriod(td))
         !  WRITE(IO_STDOUT,*) 'East TST boundary tidal number: ',  tidalnumber
        enddo
        enddo
      endif

      ! north OB
      if ( North .AND. have_g_js) then
         Vtiden(:,:)=0.
         Utiden(:,:)=0.
         Etiden(:,:)=0.
        do j = 1,3
        do td=1,tidalnumber
          Vtiden(1:ie,j)=Vtiden(1:ie,j)+                                  &
     &      OBTSTNVam(p_ioff+1:p_ioff+ie,td,j)*amsue(1:ie,j,1)*             &
     &      COS(2.0 * PI *(myTime-OBTSTNVph(p_ioff+1:p_ioff+ie,td,j))/tidalPeriod(td))

          Utiden(I_start:ie,j)=Utiden(I_start:ie,j)+                      &
     &      OBTSTNUam(p_ioff+I_start:p_ioff+ie,td,j)*amsuo(I_start:ie,j,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTNUph(p_ioff+I_start:p_ioff+ie,td,j))/tidalPeriod(td))

          Etiden(1:ie,j)=Etiden(1:ie,j)+                      &
     &      OBTSTNZam(p_ioff+1:p_ioff+ie,td,j)*weto(1:ie,j,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTNZph(p_ioff+1:p_ioff+ie,td,j))/tidalPeriod(td))
        !  WRITE(IO_STDOUT,*) 'North TST boundary tidal number: ',  tidalnumber
        enddo
        enddo
         ! WRITE(IO_STDOUT,*) 'North TST boundary Vtiden(:,1) : ',  Vtiden(1:5,1)
         ! WRITE(IO_STDOUT,*) 'North TST boundary Utiden(:,1) : ',  Utiden(1:5,1)       
         ! WRITE(IO_STDOUT,*) 'North TST boundary Etiden(:,1) : ',  Etidew(1:5,1)   
      endif

      ! south OB
      if ( South .AND. have_g_je) then
         Vtides(:,:)=0.
         Utides(:,:)=0.
         Etides(:,:)=0.
        do j = je-3,je1
        do td=1,tidalnumber
          Vtides(1:ie,j)=Vtides(1:ie,j)+                                  &
     &      OBTSTSVam(p_ioff+1:p_ioff+ie,td,j-je+4)*amsue(1:ie,j,1)*             &
     &      COS(2.0 * PI *(myTime-OBTSTSVph(p_ioff+1:p_ioff+ie,td,j-je+4))/tidalPeriod(td))

          Utides(I_start:ie,j+1)=Utides(I_start:ie,j+1)+                      &
     &      OBTSTSUam(p_ioff+I_start:p_ioff+ie,td,j-je+4)*amsuo(I_start:ie,j+1,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTSUph(p_ioff+I_start:p_ioff+ie,td,j-je+4))/tidalPeriod(td))

          Etides(1:ie,j+1)=Etides(1:ie,j+1)+                      &
     &      OBTSTSZam(p_ioff+1:p_ioff+ie,td,j-je+4)*weto(1:ie,j+1,1)* &
     &      COS(2.0 * PI * (myTime-OBTSTSZph(p_ioff+1:p_ioff+ie,td,j-je+4))/tidalPeriod(td))
        !  WRITE(IO_STDOUT,*) 'South TST boundary tidal number: ',  tidalnumber
        enddo
        enddo
      endif
#endif  /*OBCS_TST_NEST*/
      END SUBROUTINE FORC_TST_TIDES

      SUBROUTINE OBCS_APPLY_WNO
#ifdef NON_HYDROSTATIC
#ifdef OBCS_W_FORC
        IF(p_pe==p_io) THEN
          DO K=1,KEP
            if ( West ) then
               READ(IO_IN_OBCWNHW) OBWW(:,K)
            ENDIF

            if ( East ) then
               READ(IO_IN_OBCWNHE) OBEW(:,K)
            ENDIF

            if ( North ) then
               READ(IO_IN_OBCWNHN) OBNW(:,K)
            ENDIF

            if ( South ) then
               READ(IO_IN_OBCWNHS) OBSW(:,K)
            ENDIF
          ENDDO
        ENDIF

          if ( West ) then
             CALL p_bcast(OBWW,p_io)
          endif

          if ( East ) then
             CALL p_bcast(OBEW,p_io)
          endif

          if ( North ) then
             CALL p_bcast(OBNW,p_io)
          endif

          if ( South ) then
             CALL p_bcast(OBSW,p_io)
          endif

#endif   /*OBCS_W_FORC*/

          if ( West .AND. have_g_is) then
          WNO(1,:,1:ke)=OBWW(p_joff+1:p_joff+je,1:ke)*weto(1,:,:)
        endif

        if ( East .AND. have_g_ie) then
          WNO(ie,:,1:ke)=OBEW(p_joff+1:p_joff+je,1:ke)*weto(ie,:,:)
        endif

        if ( North .AND. have_g_js) then
           WNO(:,1,1:ke)=OBNW(p_ioff+1:p_ioff+ie,1:ke)*weto(:,1,:)
        endif

        if ( South .AND. have_g_je) then
           WNO(:,je,1:ke)=OBSW(p_ioff+1:p_ioff+ie,1:ke)*weto(:,je,:)
        endif
#endif /*NON_HYDROSTATIC*/

      END SUBROUTINE OBCS_APPLY_WNO

      SUBROUTINE OBCS_APPLY_WO

        if ( West .AND. have_g_is) then
          WO(1,:,1:ke)=OBWW(p_joff+1:p_joff+je,1:ke)*weto(1,:,:)
        endif

        if ( East .AND. have_g_ie) then
          WO(ie,:,1:ke)=OBEW(p_joff+1:p_joff+je,1:ke)*weto(ie,:,:)
        endif

        if ( North .AND. have_g_js) then
           WO(:,1,1:ke)=OBNW(p_ioff+1:p_ioff+ie,1:ke)*weto(:,1,:)
        endif

        if ( South .AND. have_g_je) then
           WO(:,je,1:ke)=OBSW(p_ioff+1:p_ioff+ie,1:ke)*weto(:,je,:)
        endif

      END SUBROUTINE OBCS_APPLY_WO


      SUBROUTINE OBCS_APPLY_UV
      ! west OB
      if ( West ) then
      if (have_g_is) then
        uoo(0,:,:) =obc_uw(:,:)
        voe(1,:,:) = obc_vw(:,:)
      endif
      endif

      ! east OB
      IF ( East ) THEN
      if (have_g_ie) then
        uoo(ie,:,:) = obc_ue(:,:)
        voe(ie,:,:) = obc_ve(:,:)
      endif
      endif

      ! north OB
      IF ( North ) THEN
      if (have_g_js) then
        uoo(:,1,:) = obc_un(:,:)
        voe(:,0,:)  = obc_vn(:,:)
      endif
      endif

      ! south OB
      IF ( South ) THEN
      if (have_g_je) then
        uoo(:,je,:) = obc_us(:,:)
        voe(:,je,:) = obc_vs(:,:)
      endif
      endif

      END SUBROUTINE OBCS_APPLY_UV

      SUBROUTINE OBCS_APPLY_TS
       ! west OB
      if ( West ) then
      if (have_g_is) then
        tho(1,:,:) = obc_tw(:,:)
        sao(1,:,:) = obc_sw(:,:)
      endif
      endif

      ! east OB
      IF ( East ) THEN
      if (have_g_ie) then
        tho(ie,:,:)  = obc_te(:,:)
        sao(ie,:,:) = obc_se(:,:)
      endif
      endif

      ! north OB
      IF ( North ) THEN
      if (have_g_js) then
        tho(:,1,:) = obc_tn(:,:)
        sao(:,1,:) = obc_sn(:,:)
      endif
      endif

      ! south OB
      IF ( South ) THEN
      if (have_g_je) then
        tho(:,je,:) = obc_ts(:,:)
        sao(:,je,:) = obc_ss(:,:)
      endif
      endif
      END SUBROUTINE OBCS_APPLY_TS

     SUBROUTINE UPGRADE_FLA_BT_UV
! GJQ 21.03.15 Add for TST-OBC
! Barotropic Tidal Current
#ifdef OBCS_UV_FLATHER
      ! West OB
      if ( West .AND. have_g_is) then
         USOT(I_start,:) = obc_uzw(:)*amsuo(0,:,1)
         VSET(1,:) =obc_vzw(:)*amsue(1,:,1)
      ! OB setting for FLA
       !  USO(I_start,:) = USOT(I_start,:)
       !  VSE(1,:) = VSET(1,:)
       !   WRITE(IO_STDOUT,*) 'West FLA  boundary USOT(I_start,:) : ',    USOT(I_start,2:5) 
       !   WRITE(IO_STDOUT,*) 'West FLA  boundary USO(I_start,:) : ',    USO(I_start,2:5)          
      endif

      ! East OB
      if ( East .AND. have_g_ie ) then
         USOT(ie,:) = obc_uze(:)*amsuo(ie,:,1)
         VSET(ie,:) = obc_vze(:)*amsue(ie,:,1)
      ! OB setting for FLA
        ! USO(ie,:) = USOT(ie,:)
        ! VSE(ie,:) = VSET(ie,:)
      ENDIF

      ! North OB
      if ( North .AND. have_g_js) then
        USOT(:,1) = obc_uzn(:)*amsuo(:,1,1)
        VSET(:,J_start) = obc_vzn(:)*amsue(:,0,1)
      ! OB setting for FLA
        ! USO(:,1) = USOT(:,1)
        ! VSE(:,J_start) = VSET(:,J_start)
        !  WRITE(IO_STDOUT,*) 'West FLA  boundary VSET(2:5J_start) : ',   VSET(2:5,J_start)
        !  WRITE(IO_STDOUT,*) 'West FLA  boundary VSE(2:5,J_start) : ',    VSE(2:5,J_start)
      endif

      ! South OB
      if ( South .AND. have_g_je) then
        USOT(:,je) = obc_uzs(:)*amsuo(:,je,1)
        VSET(:,je) = obc_vzs(:)*amsue(:,je,1)
      ! OB setting for FLA
        ! USO(:,je) = USOT(:,je)
        ! VSE(:,je)  = VSET(:,je)
      endif
#endif
      END SUBROUTINE UPGRADE_FLA_BT_UV


      SUBROUTINE UPGRADE_FLA_BT_ZO
#ifdef OBCS_UV_FLATHER
 ! Barotropic Tidal SeaLevel
      ! West OB
      if ( West .AND. have_g_is) then
         zetaT(1,:) = obc_zwt(:)
      endif

      ! East OB
      if ( East .AND. have_g_ie ) then
 		    zetaT(ie,:) = obc_zet(:)
      ENDIF

      ! North OB
      if ( North .AND. have_g_js) then
         zetaT(:,1) = obc_znt(:)
      endif

      ! South OB
      if ( South .AND. have_g_je) then
        zetaT(:,je) = obc_zst(:)
      endif
#endif
      END SUBROUTINE UPGRADE_FLA_BT_ZO


      SUBROUTINE UPGRADE_FLA_SUBT_ZO
#ifdef OBC_SUBTIDE
#ifdef OBCS_UV_FLATHER
 ! Barotropic SubTidal SeaLevel
      ! West OB
      if ( West .AND. have_g_is) then
         zeta_subT(1,:) = obc_zws(:)
      endif

      ! East OB
      if ( East .AND. have_g_ie ) then
 		    zeta_subT(ie,:) = obc_zes(:)
      ENDIF

      ! North OB
      if ( North .AND. have_g_js) then
         zeta_subT(:,1) = obc_zns(:)
      endif

      ! South OB
      if ( South .AND. have_g_je) then
        zeta_subT(:,je) = obc_zss(:)
      endif
#endif
#endif
      END SUBROUTINE UPGRADE_FLA_SUBT_ZO


      SUBROUTINE UPGRADE_ORL_SUBT_UV
#if defined OBC_SUBTIDE && defined ORLANSKI

        if ( West .AND. have_g_is) then
            obcf_uw = obc_usw
            obcf_vw = obc_vsw
        endif

        if ( East .AND. have_g_ie) then
              obcf_ue = obc_use
              obcf_ve = obc_vse
        endif

        if ( North .AND. have_g_js) then
               obcf_un = obc_usn
               obcf_vn = obc_vsn
        endif

        if ( South .AND. have_g_je) then
               obcf_us = obc_uss
               obcf_vs = obc_vss
        endif


#endif
      END SUBROUTINE UPGRADE_ORL_SUBT_UV


      SUBROUTINE SET_OBCS_UV(UUU0,VVV0)
#ifdef OBCS_UV_FLATHER
      REAL UUU0(1:IE-I_start+1,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,1:JE-J_start+1,KE),VVV(IE,J_start:JE,KE)
      
      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)
      
      ! The rouinte should be called before ocmodmom
        if ( West .AND. have_g_is) then
        	do k = 1,ke
        		do j = 1,je
        	  UUU(I_start,j,k) =(USOT(I_start,j) + obcf_uw(j,k))*amsuo(I_start,j,k)
        	  VVV(1,j,k) =(VSET(1,j) + obcf_vw(j,k))*amsue(1,j,k)
        	 enddo
          enddo
        endif

        if ( East .AND. have_g_ie) then
        	do k = 1,ke
        		do j = 1, je
        	  UUU(ie,j,k) =(USOT(ie,j) + obcf_ue(j,k))*amsuo(ie,j,k)
        	  VVV(ie,j,k) =(VSET(ie,j) + obcf_ve(j,k))*amsue(ie,j,k)
        	 enddo
          enddo
        endif

        if ( North .AND. have_g_js) then
        	do k = 1,ke
        		do i = 1,ie
        		 UUU(i,1,k) = (USOT(i,1) + obcf_un(i,k))*amsuo(i,1,k)
        	   VVV(i,J_start,k) = (VSET(i,J_start) + obcf_vn(i,k) )*amsue(i,J_start,k)
        	  enddo
        	  ! WRITE(IO_STDOUT,*) 'North SET_OBCS_UV k =', k
        	  ! WRITE(IO_STDOUT,*) 'VVV', VVV(2:5,J_start,k)
        	  ! WRITE(IO_STDOUT,*) 'VSET', VSET(2:5,J_start)
        	  ! WRITE(IO_STDOUT,*) 'obcf_vn', obcf_vn(2:5,k)
          enddo
        endif

        if ( South .AND. have_g_je) then
        	do k = 1,ke
        		do i = 1,ie
        	   UUU(i,je,k) =(USOT(i,je) + obcf_us(i,k))*amsuo(i,je,k)
        	   VVV(i,je,k) =(VSET(i,je) + obcf_vs(i,k))*amsue(i,je,k)
        	   
        	 !  UUU(i,je2:je1,k) = USOT(i,je2:je1)*amsuo(i,je2:je1,k)
        	 !  VVV(i,je2:je1,k) =  VSET(i,je2:je1) *amsue(i,je2:je1,k)
        	  enddo
          enddo
        endif
        
       CALL bounds_exch(UUU)
       CALL bounds_exch(VVV)

       UUU0(:,:,:)=UUU(:,:,:)
       VVV0(:,:,:)=VVV(:,:,:)
#endif      
      RETURN
      END SUBROUTINE SET_OBCS_UV


      SUBROUTINE SET_STL(ZZZ0)
     
      REAL ZZZ0(IE,JE),ZZZ(IE,JE)
      
			ZZZ(:,:) = ZZZ0(:,:)
#ifdef OBCS_TST_NEST
      ! The rouinte should be called before ocmodmom
        if ( West .AND. have_g_is) then
        	DO J=1,JE
           ZZZ(1:3,J) = ZetaT(1:3,J) * weto(1:3,J,1)
         ENDDO
         write(IO_STDOUT,*) 'SET_STL West'
         WRITE(IO_STDOUT,*)  'ZetaT(1,:)', ZetaT(1,2:10)
         WRITE(IO_STDOUT,*)  'ZZZ(1,:)', ZZZ(1,2:10)         
         WRITE(IO_STDOUT,*)  'ZetaT(2,:)', ZetaT(2,2:10)
         WRITE(IO_STDOUT,*)  'ZZZ(2,:)', ZZZ(2,2:10) 
          write(IO_STDOUT,*) 'SET_STL West'      
        endif
        
        if ( East .AND. have_g_ie) then
        	DO J=1,JE
        	ZZZ(IE2:IE,J) = ZetaT(IE2:IE,J) * weto(IE2:IE,J,1)
          ENDDO
        endif
        
        if ( North .AND. have_g_js) then
        	DO I=1,IE
        	 ZZZ(I,1:3) = ZetaT(I,1:3)* weto(I,1:3,1)
         ENDDO
        endif
        
        if ( South .AND. have_g_je) then
        	DO I=1,IE
        	  ZZZ(I,JE2:JE) = ZetaT(I,JE2:JE)* weto(I,JE2:JE,1)
         ENDDO
        endif
        
        ZZZ0(:,:) = ZZZ(:,:)
#endif
      RETURN
      END SUBROUTINE SET_STL

      SUBROUTINE PRE_TIDES
!Tidal and Subtidal signal decompostion on open boundary
!As for the open boundary, the tidal signals for rows or lines of OB, 1, 2, and 3 are required
! Subtide extraction
      ! ----1---------West Boundary-------------je----
      !      u    u    u   u    u    u   u    u    u    u      |    0  OB  (OB values setting tidal & subtidal inputs)
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    1  (OB values) -----> Updated Velocity, also needs tidal input
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    2  tidal input
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    3  tidal input
      ! N--------------inner area----------------------S |
      ! N--------------inner area----------------------S |
      !      u    u    u   u    u    u   u    u    u    u      |   ie-3  tidal input
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    ie2  tidal input
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    ie1  (OB values) -----> Updated Velocity, also needs tidal input
      !  v  p v p v p v p v p v p v p v p v p v p v   |
      !      u    u    u   u    u    u   u    u    u    u      |    ie OB values setting tidal & subtidal inputs)
      ! ----1---------East Boundary-------------je----
#ifdef OBCS_TST_NEST
! normal barotropic velocity
      if (WEST .and. have_g_is) then
      	obc_ubarw(:) = 0.
      	do j = 1,je
      		 if  (AMSUO(I_start,j,1) .eq. 1) then
      	    do k = 1,ke
      	  	   if  (AMSUO(I_start,j,k) .eq. 1) then
      		     obc_ubarw(j) = obc_ubarw(j) + obcf_uw(j,k)*DDUO(I_start,j,k)
      		     endif
      	    enddo
      	  obc_ubarw(j) = obc_ubarw(j)/DEUTO(I_start,j)
      	  endif
        enddo
      endif

      if (EAST .and. have_g_ie) then
      	obc_ubare(:) = 0.
      	do j = 1,je
      		 if  (AMSUO(ie,j,1) .eq. 1) then
      	    do k = 1,ke
      	  	   if  (AMSUO(ie,j,k) .eq. 1) then
      		     obc_ubare(j) = obc_ubare(j) + obcf_ue(j,k)*DDUO(ie,j,k)
      		     endif
      	    enddo
      	  obc_ubare(j) = obc_ubare(j)/DEUTO(ie,j)
      	  endif
        enddo
     endif

     ! ----ie---------North Boundary-------------1----
      !      v    v    v     v    v    v    v    v     v    v      |   0  OB  (OB values setting tidal & subtidal inputs)
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |   1  (OB values) -----> Updated Velocity , also needs tidal input
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |   2  tidal input
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |   3  tidal input
      ! E----------------inner area---------------------W |
      ! E----------------inner area---------------------W |
      !      v    v    v     v    v    v    v    v     v    v      |    je-3 tidal input
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |    je2  tidal input
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |    je1  (OB values) -----> Updated Velocity , also needs tidal input
      !  u  p u p u p u p u p u p u p u p u p u p u   |
      !      v    v    v    v     v    v    v    v     v    v      |    je OB  (OB values setting tidal & subtidal inputs)
      ! ----ie---------South Boundary-------------1----

      if (North .and. have_g_js) then
      	obc_vbarn(:) = 0.
      	do i = 1,ie
      		 if  (AMSUE(i,J_start,1) .eq. 1) then
      	    do k = 1,ke
      	  	   if  (AMSUE(i,J_start,k) .eq. 1) then
      		     obc_vbarn(i) = obc_vbarn(i) + obcf_vn(i,k)*DDUE(i,J_start,k)
      		     endif
      	    enddo
      	  obc_vbarn(i) = obc_vbarn(i)/DEUTE(i,J_start)
      	  endif
        enddo
      endif

     if (South .and. have_g_je) then
      	obc_vbars(:) = 0.
      	do i = 1,ie
      		 if  (AMSUE(i,je,1) .eq. 1) then
      	    do k = 1,ke
      	  	   if  (AMSUE(i,je,k) .eq. 1) then
      		     obc_vbars(i) = obc_vbars(i) + obcf_vs(i,k)*DDUE(i,je,k)
      		     endif
      	    enddo
      	  obc_vbars(i) = obc_vbars(i)/DEUTE(i,je)
      	  endif
        enddo
      endif

      if (WEST .and. have_g_is) then
         USOT(1:3,:) = Utidew(1:3,:)
         VSET(1:3,:) = Vtidew(1:3,:)
         ZetaT(1:3,:) = Etidew(1:3,:)
        ! OB normal velocity values setting
       !  USO(I_start,:) = USOT(I_start,:) + obc_ubarw(:)
       ! WRITE(IO_STDOUT,*)  'SET_TIDES WEST:   USO(I_start,:) ',   USO(I_start,2:5) 
       ! WRITE(IO_STDOUT,*)  'SET_TIDES WEST:  USOT(I_start,:)  ', USOT(I_start,2:5) 
      endif

      if (EAST .and. have_g_ie) then
         USOT(ie-3:ie1,:) = Utidee(ie-3:ie1,:)
         VSET(ie2:ie,:) = Vtidee(ie2:ie,:)
         ZetaT(ie2:ie,:) = Etidee(ie2:ie,:)
        ! OB normal velocity values setting
        ! USO(ie,:) = USOT(ie,:) + obc_ubare(:)
        ! WRITE(IO_STDOUT,*)  'EAST:  USO(ie,:) ',  USO(ie,:) 
        ! WRITE(IO_STDOUT,*)  'EAST: USOT(ie,:) ', USOT(ie,:)
      endif

      if (North .and. have_g_js) then
         VSET(:,1:3) = Vtiden(:,1:3)
         USOT(:,1:3)  = Utiden(:,1:3)
         ZetaT(:,1:3) = Etiden(:,1:3)
        ! OB normal velocity values setting
        ! VSE(:,J_start) =VSET(:,J_start)  + obc_vbarn(:)
        ! WRITE(IO_STDOUT,*)  'NORTH:  VSE(:,J_start)',  VSE(:,J_start)
        ! WRITE(IO_STDOUT,*)  'NORTH: VSET(:,J_start) ', VSET(:,J_start) 
        ! WRITE(IO_STDOUT,*)  'SET_TIDES NORTH:   VSE(2:5,J_start) ',   VSE(2:5,J_start)
        ! WRITE(IO_STDOUT,*)  'SET_TIDES NORTH:  VSET(2:5,J_start)  ', VSET(2:5,J_start)
      endif

      if (South .and. have_g_je) then
        VSET(:,je-3:je1) = Vtides(:,je-3:je1)
         USOT(:,je2:je)  = Utides(:,je2:je)
         ZetaT(:,je2:je) = Etides(:,je2:je)
        ! OB normal velocity values setting
        ! VSE(:,je) =VSET(:,je)  + obc_vbars(:)
        ! WRITE(IO_STDOUT,*)  'SOUTH: VSE(:,je)', VSE(:,je)
        ! WRITE(IO_STDOUT,*)  'SOUTH: VSET(:,je)', VSET(:,je)
      endif
#endif  /*OBCS_TST_NEST*/
     END SUBROUTINE PRE_TIDES


     SUBROUTINE OBCS_FLATHER
#ifdef OBCS_UV_FLATHER
      REAL bry_val
! ------------------ U_bt --------------------------------
     ! Western Edge -- U_bt
       if (WEST .and. have_g_is) then
            if (North .and. have_g_js) then
        	     j_begin = 2    ! NorthWest Corner Excluded
        	     j_end = je
           else if (South .and. have_g_je) then
           	   j_begin = 1
           	   j_end = je1   ! SouthWest Corner Excluded
           	else
               j_begin = 1
           	   j_end = je
           endif
           do j=j_end,j_begin,-1
             if (AMSUO(1,J,1) .eq. 1) THEN
               bry_val = USOT(0,j)   + obc_ubarw(j)   ! BT_tide + BT_subtide
               ! cff = 1.0 / (DEUTO(0,j) + 0.5*(ZOOLD(1,j)+ ZOOLD(2,j)))  ! ATTEN: We use DEUTO!
               cff = 1.0 / (0.5*(DEPTO(1,j)+ZOOLD(1,j) + DEPTO(2,j)+ ZOOLD(2,j) ) )

               Cx = SQRT(G*cff)
               USO(1,j) = bry_val - Cx * (0.5*(ZOOLD(1,j)+ ZOOLD(2,j))-(zeta_subT(1,j)+zetaT(1,j)))
               USO(1,j) = USO(1,j) * AMSUO(1,j,1)
             endif
          end do
       end if

     ! Eastern Edge -- U_bt
       if (EAST .and. have_g_ie) then
            if (North .and. have_g_js) then
        	     j_begin = 2  ! NorthEast Corner Excluded
        	     j_end = je
           else if (South .and. have_g_je) then
           	   j_begin = 1
           	   j_end = je1  ! SouthEast Corner Excluded
           	else
               j_begin = 1
           	   j_end = je
           endif

           do j=j_end,j_begin,-1
               if (AMSUO(IE1,J,1) .eq. 1) THEN
                 bry_val =  USOT(ie,j) + obc_ubare(j)
                 ! cff = 1.0 / (DEUTO(ie,j) + 0.5*(ZOOLD(ie,j)+ZOOLD(ie1,j)))
                 cff = 1.0 / (0.5*(DEPTO(ie,j)+ZOOLD(ie,j) + DEPTO(ie1,j)+ ZOOLD(ie1,j) ) )

                 Cx = SQRT(G*cff)
                 USO(ie1,j) = bry_val + Cx * (0.5*(ZOOLD(ie,j)+ZOOLD(ie1,j))-(zeta_subT(ie,j)+zetaT(ie,j)))
                 USO(ie1,j) = USO(ie1,j) * AMSUO(ie1,j,1)
               end if
        end do
       end if


! Southern Edge -- U_bt *ATTEN: This is Chapman BC.
        if (ChpmSouth) then
        if (SOUTH .and. have_g_je )  then
            if (WEST .and. have_g_is) then
        	     i_begin = 2 ! SouthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1
           	   i_end = ie2  ! SouthEast Corner Excluded
           	else
               i_begin = 1
           	   i_end = ie1
           endif

           do i=i_begin,i_end
              if (AMSUO(I,JE,1) .eq. 1) THEN
                 cff  = dpyo(i,je1)  ! DPYO = DT / DLYU
              ! ATTEN: JE1 is for NORTH - SOUTH.
              ! C = - SQRT(G*D)
      !          cff1 = SQRT(G*(DEUTO(i,je1) +                           &
      ! &               0.5 * (ZOOLD(i,je1) + ZOOLD(i+1,je1))))
                 cff1 = SQRT(G*0.5*(DEPTO(i,je1) +ZOOLD(i,je1)     &
     &                     +DEPTO(i+1,je1)+ ZOOLD(i+1,je1)))
                 Ce   = cff * cff1
                 cff2 = 1.0 / (1.0 + Ce)
                 USO(i,je) = cff2 * (USOOLD(i,je) +                          &
     &                         Ce*USO(i,je1))
                 USO(i,je) = USO(i,je) * AMSUO(i,je,1)
              end if
          end do

         endif
        endif


! Northern Edge -- U_bt
       if (ChpmNorth) then
       if (NORTH .and. have_g_js) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2  ! NorthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1
           	   i_end = ie2   ! NorthEast Corner Excluded
           	else
               i_begin = 1
           	   i_end = ie1
           endif
           do i=i_begin,i_end
              if (AMSUO(I,1,1) .eq. 1) THEN
                cff = dpyo(i,2) ! DPYO = DT / DLYU
                ! ATTEN: 2 is for NORTH - SOUTH.
                ! C = SQRT(G*D)
!               cff1 = SQRT(G*(DEUTO(i,2) +                               &
!    &               0.5 * (ZOOLD(i,2) + ZOOLD(i+1,2))))
                cff1 = SQRT(G*0.5*(DEPTO(i,2) +ZOOLD(i,2)     &
     &                     +DEPTO(i+1,2)+ ZOOLD(i+1,2)))
                Ce   = cff * cff1
                cff2 = 1.0 / (1.0 + Ce)
                USO(i,1) = cff2 * (USOOLD(i,1) +                           &
     &                         Ce*USO(i,2))
                USO(i,1) = USO(i,1) * AMSUO(i,1,1)
            end if
         end do
        endif
       endif

! ------------------ V_bt --------------------------------

! Southern Edge -- V_bt
      if (SOUTH .and. have_g_je) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2 ! SouthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1  ! SouthEast Corner Excluded
           	   i_end = ie1
           	else
               i_begin = 1
           	   i_end = ie
           endif
          do i=i_begin,i_end
              if (AMSUE(I,JE1,1) .eq. 1) THEN
               bry_val = obc_vbars(i) + VSET(i,je)
!              cff = 1.0 / ( DEUTE(i,je) + 0.5*(ZOOLD(i,je) +ZOOLD(i,je1)) )
               cff = 1.0 /(0.5* ( DEPTO(i,je) + ZOOLD(i,je) +DEPTO(i,je1) + ZOOLD(i,je1) ) )

               Ce  = SQRT(G*cff)
               VSE(i,je1) = bry_val - Ce * (0.5*(ZOOLD(i,je) +ZOOLD(i,je1)) - (zeta_subT(i,je)+zetaT(i,je)))
               VSE(i,je1) = VSE(i,je1) * AMSUE(i,je1,1)
              end if
          end do
      end if

! Northern Edge -- V_bt
      if (NORTH .and. have_g_js) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2   ! SouthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1
           	   i_end = ie1   ! SouthEast Corner Excluded
           	else
               i_begin = 1
           	   i_end = ie
           endif
           do i=i_begin,i_end
               if (AMSUE(I,1,1) .eq. 1) THEN
                 bry_val = obc_vbarn(i) + VSET(i,0)
                 ! cff = 1.0 / ( DEUTE(i,0) + 0.5*(ZOOLD(i,1) + ZOOLD(i,2)) )
                 cff = 1.0 /(0.5* ( DEPTO(i,1) + ZOOLD(i,1) +DEPTO(i,2) + ZOOLD(i,2) ) )

                 Ce  = SQRT(G*cff)
                 VSE(i,1) = bry_val + Ce * (0.5*(ZOOLD(i,1) + ZOOLD(i,2)) - (zeta_subT(i,0)+zetaT(i,0)))
                 VSE(i,1) = VSE(i,1) * AMSUE(i,1,1)
               end if
           end do
      end if


! Western Edge -- V_bt
      if (ChpmWest) then
      if (WEST .and. have_g_is) then
      	if (North .and. have_g_js) then
      		  j_begin = 2  ! NorthWest Corner Excluded
      		  j_end = je
      	else if (South .and. have_g_je) then
      		  j_begin = 1
      		  j_end = je2  ! SouthWest Corner Excluded
        else
      		  j_begin = 1
      		  j_end = je1
      	endif

         do j=j_end,j_begin,-1
           if (AMSUE(1,J,1) .eq. 1) THEN
           cff  = DTDXUE(2,j) ! DTXDUE = DT / DLXV
           ! ATTEN: 2 is for NORTH - SOUTH.
           ! C = - SQRT(G*D)
 !          cff1 = SQRT(G*(DEUTE(2,j) +                               &
 !    &               0.5 * (ZOOLD(2,j ) + ZOOLD(2,j+1))))
           cff1 = SQRT(G*0.5*(DEPTO(2,j) +ZOOLD(2,j)     &
     &                     +DEPTO(2,j+1)+ ZOOLD(2,j+1)))

           Cx   = cff * cff1
           cff2 = 1.0 / (1.0 + Cx)
           VSE(1,j) = cff2 * (VSEOLD(1,j) +                             &
     &                         Cx*VSE(2,j))
           VSE(1,j) = VSE(1,j) * AMSUE(1,j,1)
           end if
        end do
      end if
      endif

! Eastern Edge -- V_bt
      if (ChpmEast) then
      if (EAST .and. have_g_ie) then
      	if (North .and. have_g_js) then
      		  j_begin = 2   ! NorthEast Corner Excluded
      		  j_end = je
      	else if (South .and. have_g_je) then
      		  j_begin = 1
      		  j_end = je2    ! SouthEast Corner Excluded
        else
      		  j_begin = 1
      		  j_end = je1
      	endif

         do j=j_end,j_begin,-1
           if (AMSUE(IE,J,1) .eq. 1) THEN
           cff  = DTDXUE(ie1,j) ! DTXDUE = DT / DLXV
           ! ATTEN: JE1 is for NORTH - SOUTH.
           ! C =  SQRT(G*D)
!          cff1 = SQRT(G*(DEUTE(ie1,j) +                             &
!     &               0.5 * (ZOOLD(ie1,j ) + ZOOLD(ie1,j+1))))
           cff1 = SQRT(G*0.5*(DEPTO(ie1,j) +ZOOLD(ie1,j)     &
     &                     +DEPTO(ie1,j+1)+ ZOOLD(ie1,j+1)))

           Cx   = cff * cff1
           cff2 = 1.0 / (1.0 + Cx)
           VSE(ie,j) = cff2 * (VSEOLD(ie,j) +                           &
     &                         Cx*VSE(ie1,j))
           VSE(ie,j) = VSE(ie,j) * AMSUE(ie,j,1)
           end if
         end do
      end if
      endif
      
      if (Cornersmean) then
! ------------------ BOUNDARY CORNERS BEGIN --------------------------------
! NorthWest Corner
      if (have_g_is .and. have_g_js) then
         if (WEST .and. NORTH) then
         	  if (WETO(1,1,1) .eq. 1) then
              USO(I_start+1,1) = 0.5*(USO(I_start+2,1) +USO(I_start+1,2))*AMSUO(I_start+1,1,1)
              VSE(1,1) = 0.5*(VSE(1,2) +VSE(2,1) )*AMSUE(1,1,1)
           endif
        endif
     endif
! SouthWest Corner
      if (have_g_is .and. have_g_je) then
      	 if (WEST .and. SOUTH) then
      	 	 if (WETO(1,je,1) .eq. 1) then
               USO(I_start+1,je) = 0.5*(USO(I_start+2,je)+USO(I_start+1,je1) )*AMSUO(I_start+1,je,1)
               VSE(1,je1) = 0.5*(VSE(2,je1) +VSE(1,je2) )*AMSUE(1,je1,1)
           endif
        endif
     endif
! NorthEast Corner
      if (have_g_ie .and. have_g_js) then
         if (EAST .and. NORTH) then
         	  if (WETO(ie,1,1) .eq. 1) then
               USO(ie1,1) = 0.5*(USO(ie2,1) +USO(ie1,2) )*AMSUO(ie1,1,1)
               VSE(ie,1) = 0.5*(VSE(ie1,1) +VSE(ie,2) )*AMSUE(ie,1,1)
           endif
        endif
     endif
! SouthEast Corner
      if (have_g_ie .and. have_g_je) then
         if (EAST .and. SOUTH) then
         	  if (WETO(ie,je,1) .eq. 1) then
              USO(ie1,je) = 0.5*(USO(ie2,je) +USO(ie1,je1) )*AMSUO(ie1,je,1)
              VSE(ie,je1) = 0.5*(VSE(ie,je2) +VSE(ie1,je1) )*AMSUE(ie,je1,1)
          endif
        endif
     endif
! ------------------ BOUNDARY CORNERS END --------------------------------
    endif
    
    
      CALL bounds_exch(USO)
      CALL bounds_exch(VSE)


#endif /*OBCS_UV_FLATHER*/

      END SUBROUTINE OBCS_FLATHER


      SUBROUTINE OBCS_CHAPMAN
#ifdef OBCS_UV_FLATHER
!GJQ 21.03.22 WHY there's not a i_step = 0 if...
!             Because the ZOOLD can be start from ZERO.
      ZONEW = 0.0
!H.H 01.05.23  Implicit Chapman boundary condition
!Refered from Zeta Updated method in zetabc.F of ROMS

! Western Edge
      if (WEST .and. have_g_is) then
            if (North .and. have_g_js) then
        	     j_begin = 2    ! NorthWest Corner Excluded
        	     j_end = je
           else if (South .and. have_g_je) then
           	   j_begin = 1
           	   j_end = je1   ! SouthWest Corner Excluded
           	else
               j_begin = 1
           	   j_end = je
           endif
         do j=j_end,j_begin,-1
           if (WETO(1,J,1) .eq. 1) THEN
            cff  = DTDXPO(2,j)
            cff1 = SQRT(G*(DEPTO(2,j) + ZOOLD(2,j)))
            Cx   = cff * cff1
            cff2 = 1.0 / (1.0 + Cx)
            ZONEW(1,j) = cff2*(ZOOLD(1,j) + Cx * (ZO(2,j) + Z1O(2,j)))
            Z1O(1,j) = ZONEW(1,j) - ZO(1,j)
            Z1O(1,j) = Z1O(1,j) * WETO(1,j,1)
           end if
         end do
      end if

! Eastern Edge

      if (EAST .and. have_g_ie) then
            if (North .and. have_g_js) then
        	     j_begin = 2  ! NorthEast Corner Excluded
        	     j_end = je
           else if (South .and. have_g_je) then
           	   j_begin = 1
           	   j_end = je1  ! SouthEast Corner Excluded
           	else
               j_begin = 1
           	   j_end = je
           endif
         do j=j_end,j_begin,-1
           if (WETO(IE,J,1) .eq. 1) THEN
            cff  = DTDXPO(ie1,j)
            cff1 = SQRT(G*(DEPTO(ie1,j) + ZOOLD(ie1,j)))
            Cx   = cff * cff1
            cff2 = 1.0 / (1.0 + Cx)
            ZONEW(ie,j) = cff2*(ZOOLD(ie,j) + Cx * (ZO(ie1,j) +         &
     &                          Z1O(ie1,j))  )
            Z1O(ie,j) = ZONEW(ie,j) - ZO(ie,j)
            Z1O(ie,j) = Z1O(ie,j) * WETO(ie,j,1)
           end if
         end do
      end if

! Southern Edge

      if (SOUTH .and. have_g_je) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2 ! SouthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1  ! SouthEast Corner Excluded
           	   i_end = ie1
           	else
               i_begin = 1
           	   i_end = ie
           endif
         do i=i_begin,i_end
           if (WETO(I,JE,1) .eq. 1) THEN
            cff  = DTDYO(i,je1)
            cff1 = SQRT(G*(DEPTO(i,je1) + ZOOLD(i,je1)))
            Ce   = cff * cff1
            cff2 = 1.0 / (1.0 + Ce)
            ZONEW(i,je) = cff2*(ZOOLD(i,je) + Ce * (ZO(i,je1) +         &
     &                          Z1O(i,je1)) )
            Z1O(i,je) = ZONEW(i,je) - ZO(i,je)
            Z1O(i,je) = Z1O(i,je) * WETO(i,je,1)
           end if
         end do
      end if

! Northern Edge
      if (NORTH .and. have_g_js) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2 ! NorthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1  ! NorthEast Corner Excluded
           	   i_end = ie1
           	else
               i_begin = 1
           	   i_end = ie
           endif
         do i=i_begin,i_end
           if (WETO(I,1,1) .eq. 1) THEN
            cff  = DTDYO(i,2)
            cff1 = SQRT(G*(DEPTO(i,2) + ZOOLD(i,2)))
            Ce   = cff * cff1
            cff2 = 1.0 / (1.0 + Ce)
            ZONEW(i,1) = cff2*(ZOOLD(i,1) + Ce * (ZO(i,2) + Z1O(i,2))  )
            Z1O(i,1) = ZONEW(i,1) - ZO(i,1)
            Z1O(i,1) = Z1O(i,1) * WETO(i,1,1)
           end if
         end do
      end if

! ------------------ BOUNDARY CORNERS --------------------------------
! NorthWest Corner
      if (have_g_is .and. have_g_js) then
         if (WEST .and. NORTH) then
         	  if (WETO(1,1,1) .eq. 1) then
               Z1O(1,1) = 0.5*(Z1O(2,1) +Z1O(1,2) )*WETO(1,1,1)
           endif
        endif
     endif
! SouthWest Corner
      if (have_g_is .and. have_g_je) then
      	 if (WEST .and. SOUTH) then
      	 	 if (WETO(1,je,1) .eq. 1) then
              Z1O(1,je) = 0.5*(Z1O(2,je) +Z1O(1,je1) )*WETO(1,je,1)
           endif
        endif
     endif
! NorthEast Corner
      if (have_g_ie .and. have_g_js) then
         if (EAST .and. NORTH) then
         	  if (WETO(ie,1,1) .eq. 1) then
               Z1O(ie,1) = 0.5*(Z1O(ie1,1) +Z1O(ie,2) )*WETO(ie,1,1)
           endif
        endif
     endif
! SouthEast Corner
      if (have_g_ie .and. have_g_je) then
         if (EAST .and. SOUTH) then
         	  if (WETO(ie,je,1) .eq. 1) then
              Z1O(ie,je) = 0.5*(Z1O(ie1,je) +Z1O(ie,je1) )*WETO(ie,je,1)
          endif
        endif
     endif
      CALL bounds_exch(Z1O)
#endif /*OBCS_UV_FLATHER*/
      END SUBROUTINE OBCS_CHAPMAN

! ---------------------- Orlanski Radiation OBCS --------------------------- !
! Notes: : Radiation: Flather (BT) + Orlanski (BC)
! ------------------------------------------------------------------------------------- !
      SUBROUTINE ORLANSKI_UV(UUU0,VVV0,i_step)
#ifdef ORLANSKI
! This is the original Radiation obc for BC current.
      REAL UUU0(1:IE-I_start+1,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,1:JE-J_start+1,KE),VVV(IE,J_start:JE,KE)
      REAL M3obc_in, M3obc_out,M3nudgcof
      INTEGER obcfac 
      ! M3obc_in = 4.0 * ( 1.0/(1.0*86400.0) ) ! ==OBC_fac*(1/OBCCOF(DAYS)*86400)
      ! OBCcof: Nudging/Relaxtion time scale.
      ! OBCfac: Factor between passive(outflow) and active(inflow) open
      !         boundary conditions, If OBCFAC > 1, nudging on inflow is
      !         stronger
      !         than on outflow (which is recommended).
      ! M3obc_out = ( 1.0/(1.0*86400.0) )

      obcfac = 4
      M3nudgcof = 1.0*86400
      M3obc_out = M3nudgcof
      M3obc_in = obcfac * M3obc_out
      
      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)


! -----------------------  U-bcs for 4 directions ---------------------------
!----------------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.(0 = I_start =
!  Istr,je1 = Jstr)
!----------------------------------------------------------------------------
 !----------------------------------------------------------------------------
  ! ATTEN: The mpiom's je is from north to south,oppsite to Roms.
  ! gradu_wb(0:2,0:je)
  ! Partial U / Partial n = 0

! Western edge, all BCs in the ROMS code.
      if (WEST .and. have_g_is) then
       gradu_wb=0.0
       ! uoo_wb0 = UUU(0:2,:,:)
       uoo_wb0(0:2,:,:) = UUU(1:3,:,:)
       do k=1,ke  ! k
         do j=je1,1,-1
            if (AMSUO(0,J,k)  .eq. 1) THEN
           ! dUdy
           ! j = je  gradu_wb = 0
           ! j = 0   gradu_wb = 0  apply !
            gradu_wb(0,j) = uoo_wb1(0,j,k) - uoo_wb1(0,j+1,k)
            gradu_wb(1,j) = uoo_wb1(1,j,k) - uoo_wb1(1,j+1,k)
            end if
         enddo
  ! Attenetion:  the corners excluded
            do j=je1,2,-1
               if (AMSUO(0,J,k)  .eq. 1) THEN
               dUdt=uoo_wb1(1,j,k) - uoo_wb0(1,j,k)
            !   dUdx=uoo_wb0(0,j,k) - uoo_wb0(1,j,k)
               dUdx=uoo_wb0(1,j,k) - uoo_wb0(2,j,k)  ! upstream shceme

               if ((dUdt*dUdx) .lt. 0.) dUdt=0.
               if ((dUdt*(gradu_wb(1,j)+gradu_wb(1,j-1))) .gt. 0.) then  ! upstream shceme
                  dUde = gradu_wb(1,j)
               else
                  dUde = gradu_wb(1,j-1)
               endif
               cff = max(dUdx*dUdx+dUde*dUde,obcerr)
               Cx = dUdt*dUdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
                obc_uw(j,k) = (cff*uoo_wb1(0,j,k)+                       &
     &                        Cx *uoo_wb0(1,j,k)-                       &
     &                        MAX(Ce,0.)*gradu_wb(0,j  )-                &
     &                        MIN(Ce,0.)*gradu_wb(0,j-1))/               &
     &                        (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dUdt*dUdx) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif


#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_uw(j,k) = (1.0 - phi)*obc_uw(j,k) +                  &
     &                       phi*(obcf_uw(j,k) + USOT(0,j))
#else
               obc_uw(j,k) = obc_uw(j,k) +                              &
     &                       tau*((obcf_uw(j,k)+USOT(0,j))-             &
     &                            uoo_wb1(0,j,k) )
#endif


#endif /*RADIATION_NUDGING*/
               endif
           enddo
       enddo  ! k
       uoo_wb1=uoo_wb0
      end if


!---------------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
! (ie1 = Iend, je1 = Jstr)
! gradu_eb(ie1:ie,0:je)
! Partial U / Partial n = 0
!---------------------------------------------------------------------------
      if (EAST .and. have_g_ie) then
         gradu_eb=0.0
        ! uoo_eb0=UUU(ie2:ie,:,:)
         uoo_eb0(ie2:ie,:,:)=UUU(ie-3:ie-1,:,:)
         do k=1,ke
            do j=je1,1,-1
               if (AMSUO(IE,J,K) .eq. 1) THEN
               ! dUdy
               ! j = je  gradu_wb = 0
               ! j = 0   gradu_wb = 0  apply !
               gradu_eb(ie1,j) = uoo_eb1(ie1,j,k) - uoo_eb1(ie1,j+1,k)
               gradu_eb(ie, j) = uoo_eb1(ie, j,k) - uoo_eb1(ie ,j+1,k)
               end if
            end do
  ! Attenetion:  the corners excluded
            do j=je1,2,-1
               if (AMSUO(IE,J,K) .eq. 1) THEN
               dUdt=uoo_eb1(ie1,j,k) - uoo_eb0(ie1,j,k)
               dUdx=uoo_eb0(ie1,j,k) - uoo_eb0(ie2,j,k) ! upstream shceme

            if ((dUdt*dUdx) .lt. 0.) dUdt=0.
            if ((dUdt*(gradu_eb(ie1,j)+gradu_eb(ie1,j-1))) .gt. 0.) then  ! upstream shceme
               dUde = gradu_eb(ie1,j)
            else
               dUde = gradu_eb(ie1,j-1)
            end if
            cff = max(dUdx*dUdx+dUde*dUde,obcerr)
            Cx = dUdt*dUdx
#ifdef RADIATION_2D
            Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
            Ce = 0.
#endif /*RADIATION_2D*/
            obc_ue(j,k) = (cff*uoo_eb1(ie, j,k)+                        &
     &                     Cx *uoo_eb0(ie1,j,k)-                        &
     &                     MAX(Ce,0.)*gradu_eb(ie,j  )-                  &
     &                     MIN(Ce,0.)*gradu_eb(ie,j-1))/                 &
     &                     (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dUdt*dUdx) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_ue(j,k) = (1.0 - phi)*obc_ue(j,k) +                  &
     &                       phi * (obcf_ue(j,k) + USOT(ie,j))
#else
               obc_ue(j,k) = obc_ue(j,k) +                              &
     &                       tau*((obcf_ue(j,k) + USOT(ie,j))-          &
     &                            uoo_eb1(ie,j,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
            end do
         end do
         uoo_eb1=uoo_eb0
      endif

!----------------------------------------------------------------------------
!  Later boundary conditions at the southern edge.(1 = I_start = IstrU,
!  ie1 = Iend, je-1 = Jstr)
!----------------------------------------------------------------------------
      if (TanSouth) then
      if (SOUTH .and. have_g_je) then
         gradu_sb=0.0
         uoo_sb0=UUU(:,je2:je,:)
         do k=1,ke
            do i=I_start,ie1
               if (AMSUO(I,JE,K) .eq. 1) THEN
               	! dUdx
               gradu_sb(i,je) =uoo_sb1(i+1,je ,k) - uoo_sb1(i,je ,k)
               gradu_sb(i,je1)=uoo_sb1(i+1,je1,k) - uoo_sb1(i,je1,k)
               end if
            end do
! Attentions: the corners excluded
               do i=I_start+2,ie2
                  if (AMSUO(I,JE,K) .eq. 1) THEN
                  dUdt=uoo_sb1(i,je1,k) - uoo_sb0(i,je1,k)
                  dUde=uoo_sb0(i,je1,k) - uoo_sb0(i,je2,k)  ! upstream shceme

                  if ((dUdt*dUde) .lt. 0.) dUdt=0.
                  if ((dUdt*(gradu_sb(i-1,je1)+gradu_sb(i,je1))) .gt. 0.)then ! upstream shceme
                     dUdx = gradu_sb(i-1,je1)
                  else
                     dUdx = gradu_sb(i  ,je1)
                  end if
                  cff=MAX(dUdx*dUdx+dUde*dUde,obcerr)
                  Ce=dUdt*dUde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dUdt*dUdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_us(i,k) = (cff*uoo_sb1(i,je, k)+                  &
     &                           Ce *uoo_sb0(i,je1,k)-                  &
     &                           MAX(Cx,0.)*gradu_sb(i-1,je)-            &
     &                           MIN(Cx,0.)*gradu_sb(i  ,je))/           &
     &                           (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dUdt*dUde) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_us(i,k) = (1.0 - phi)*obc_us(i,k) +                  &
     &                       phi * (obcf_us(i,k) + USOT(i,je))
#else
               obc_us(i,k) = obc_us(i,k) +                              &
     &                       tau*((obcf_us(i,k) + USOT(i,je))-          &
     &                            uoo_sb1(i,je,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
               end do
         end do
         uoo_sb1=uoo_sb0
      endif
      endif

!-----------------------------------------------------------------------------
!  Later boundary conditions at the northern edge.(1 = I_start = IstrU,
!  2 = Jend)
!-----------------------------------------------------------------------------
      if (TanNorth) then
      if (NORTH .and. have_g_js) then
         gradu_nb=0.0
         uoo_nb0=UUU(:,1:3,:)  ! allocate(uoo_nb0(I_start:ie,1:3,ke))
         do k=1,ke
            do i=I_start,ie1
               if (AMSUO(I,1,K) .eq. 1) THEN
               	! dUdx
               gradu_nb(i,2) = uoo_nb1(i+1,2,k) - uoo_nb1(i,2,k)
               gradu_nb(i,1) = uoo_nb1(i+1,1,k) - uoo_nb1(i,1,k)
               end if
            end do
! Attentions: the corners excluded
               do i=I_start+2,ie2
                  if (AMSUO(I,1,K) .eq. 1) THEN
                  dUdt=uoo_nb1(i,2,k) - uoo_nb0(i,2,k)
                  dUde=uoo_nb0(i,2,k) - uoo_nb0(i,3,k)  ! upstream shceme

                  if ((dUdt*dUde) .lt. 0.) dUdt=0.
                  if ((dUdt*(gradu_nb(i-1,2)+gradu_nb(i,2))) .gt. 0.) then
                     dUdx = gradu_nb(i-1,2)
                  else
                     dUdx = gradu_nb(i  ,2)
                  end if
                  cff=MAX(dUdx*dUdx+dUde*dUde,obcerr)
                  Ce=dUdt*dUde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dUdt*dUdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                 obc_un(i,k) = (cff*uoo_nb1(i,1,k)+                     &
     &                          Ce *uoo_nb0(i,2,k)-                     &
     &                          MAX(Cx,0.)*gradu_nb(i-1,1)-              &
     &                          MIN(Cx,0.)*gradu_nb(i  ,1))/             &
     &                          (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dUdt*dUde) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_un(i,k) = (1.0 - phi)*obc_un(i,k) +                  &
     &                       phi * (obcf_un(i,k) + USOT(i,1))
#else
               obc_un(i,k) = obc_un(i,k) +                              &
     &                       tau*((obcf_un(i,k)+USOT(i,1))-             &
     &                            uoo_nb1(i,1,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
               end do
         end do
         uoo_nb1=uoo_nb0
      endif
      endif


! U_bc is done! NEXT is V_bc...

!-------------------------  V-bcs for 4 directions
!---------------------------
!-----------------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.(2 = Istr,ie =
!  Iend+1, je = Jstr)
!  gradv_sb(1:ie+1,je1:je)
!-----------------------------------------------------------------------------
! Southern edge
      if (SOUTH .and. have_g_je) then
         gradv_sb=0.0
         ! voe_sb0=VVV(:,je2:je,:)
         voe_sb0(:,je2:je,:)=VVV(:,je-3:je-1,:)
         do k=1,ke
            do i=2,ie
               if (AMSUE(I,JE,K) .eq. 1) THEN
               	! dVdx
               	! i = 1 gradv_sb = 0
               	! i = ie+1  gradv_sb = 0
               gradv_sb(i,je)   = voe_sb1(i,je, k) - voe_sb1(i-1,je, k)
               gradv_sb(i,je1)  = voe_sb1(i,je1,k) - voe_sb1(i-1,je1,k)
               end if
            end do
! Attentions: the corners excluded
               do i=2,ie1
                  if (AMSUE(I,JE,K) .eq. 1) THEN
                  dVdt=voe_sb1(i,je1,k) - voe_sb0(i,je1,k)
                  dVde=voe_sb0(i,je1,k) - voe_sb0(i,je2,k)  ! upstream shceme

               if ((dVdt*dVde).lt.0.) dVdt=0.
               if ((dVdt*(gradv_sb(i,je1)+gradv_sb(i+1,je1))).gt.0.) then ! upstream shceme
                  dVdx=gradv_sb(i,je1)
               else
                  dVdx=gradv_sb(i+1,je1)
               end if
               cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
               Ce=dVdt*dVde
#ifdef RADIATION_2D
               Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#else
               Cx=0.
#endif /*RADIATION_2D*/
               obc_vs(i,k) = (cff*voe_sb1(i,je,k)+                      &
     &                        Ce *voe_sb0(i,je1,k)-                     &
     &                        MAX(Cx,0.)*gradv_sb(i,je)-                 &
     &                        MIN(Cx,0.)*gradv_sb(i+1,je))/              &
     &                        (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dVdt*dVde) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_vs(i,k) = (1.0 - phi)*obc_vs(i,k) +                  &
     &                       phi * (obcf_vs(i,k) + VSET(i,je))
#else
               obc_vs(i,k) = obc_vs(i,k) +                              &
     &                       tau*((obcf_vs(i,k) + VSET(i,je))-          &
     &                            voe_sb1(i,je,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
               end do
         end do
         voe_sb1=voe_sb0
      endif

!---------------------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.(2 = Istr, 1 = Jend,
!  0 = Jend-1)
! gradv_nb(1:ie+1,0:2)
!---------------------------------------------------------------------------------
! northern edge
      if (NORTH .and. have_g_js) then
         gradv_nb=0.0
         ! voe_nb0=VVV(:,0:2,:)
         voe_nb0(:,0:2,:)=VVV(:,1:3,:)
         do k=1,ke
            do i=2,ie
               if (AMSUE(I,0,K) .eq. 1) THEN
               	! i = 1 gradv_nb = 0
               	! i = ie+1  gradv_nb = 0
               gradv_nb(i,1) = voe_nb1(i,1,k) - voe_nb1(i-1,1,k)
               gradv_nb(i,0) = voe_nb1(i,0,k) - voe_nb1(i-1,0,k)
               end if
            end do
! Attentions: the corners excluded
                do i=2,ie1
                  if (AMSUE(I,0,K) .eq. 1) THEN
                  dVdt=voe_nb1(i,1,k) - voe_nb0(i,1,k)
                  dVde=voe_nb0(i,1,k) - voe_nb0(i,2,k)  ! upstream shceme

                  if ((dVdt*dVde).lt.0.) dVdt=0.
                  if ((dVdt*(gradv_nb(i,1)+gradv_nb(i+1,1))).gt.0.) then ! upstream shceme
                     dVdx=gradv_nb(i  ,1)
                  else
                     dVdx=gradv_nb(i+1,1)
                  end if
                  cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
                  Ce=dVdt*dVde
#ifdef RADIATION_2D
                  Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#else
                  Cx=0.0
#endif /*RADIATION_2D*/
                  obc_vn(i,k) = (cff*voe_nb1(i,0,k)+                    &
     &                           Ce *voe_nb0(i,1,k)-                    &
     &                           MAX(Cx,0.)*gradv_nb(i,  0)-             &
     &                           MIN(Cx,0.)*gradv_nb(i+1,0))/            &
     &                           (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dVdt*dVde) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_vn(i,k) = (1.0 - phi)*obc_vn(i,k) +                  &
     &                       phi * (obcf_vn(i,k) + VSET(i,0))
#else
               obc_vn(i,k) = obc_vn(i,k) +                              &
     &                       tau*((obcf_vn(i,k) + VSET(i,0))-           &
     &                            voe_nb1(i,0,k) )
#endif

#endif /*RADIATION_NUDGING*/
              end if
               end do
         end do
         voe_nb1=voe_nb0
      end if

!--------------------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.(je = JstrV-1, je1 =
!  Jend, 2 = Istr)
!  gradv_wb(1:3,1:je)
!--------------------------------------------------------------------------------
      if (TanWest) then
      if (WEST .and. have_g_is) then
         gradv_wb=0.0
         voe_wb0=VVV(1:3,:,:)
         do k=1,ke
            do j=je,J_start+1,-1 ! ATTEN: The direction is from S -> N.
               if (AMSUE(1,J,K) .eq. 1) THEN
               	! dVdy
               gradv_wb(1,j) =  voe_wb1(1,j-1,k) - voe_wb1(1,j,k)
               gradv_wb(2,j) =  voe_wb1(2,j-1,k) - voe_wb1(2,j,k)
               end if
            end do
! Attentions: the corners excluded
               do j=je2,J_start+2,-1
                  if (AMSUE(1,J,K) .eq. 1) THEN
                  dVdt=voe_wb1(2,j,k)-voe_wb0(2,j,k)
                  dVdx=voe_wb0(2,j,k)-voe_wb0(3,j,k)  ! upstream shceme

                  if ((dVdt*dVdx).lt.0.) dVdt=0.
                  if ((dVdt*(gradv_wb(2,j+1)+gradv_wb(2,j))).gt.0.) then ! upstream shceme
                     dVde=gradv_wb(2,j+1)
                  else
                     dVde=gradv_wb(2,j  )
                  end if
                     cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
                     Cx=dVdt*dVdx
#ifdef RADIATION_2D
                     Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
                     Ce=0.
#endif /*RADIATION_2D*/
                     obc_vw(j,k) = (cff*voe_wb1(1,j,k)+                 &
     &                             Cx *voe_wb0(2,j,k)-                  &
     &                             MAX(Ce,0.)*gradv_wb(1,j+1)-           &
     &                             MIN(Ce,0.)*gradv_wb(1,j  ))/          &
     &                             (cff+Cx)

#ifdef RADIATION_NUDGING
               if ((dVdt*dVdx) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_vw(j,k) = (1.0 - phi)*obc_vw(j,k) +                  &
     &                       phi * (obcf_vw(j,k) + VSET(1,j))
#else
               obc_vw(j,k) = obc_vw(j,k) +                              &
     &                       tau*((obcf_vw(j,k) + VSET(1,j))-           &
     &                            voe_wb1(1,j,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
               end do
         end do
         voe_wb1=voe_wb0
      end if
      endif

!-------------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.(ie1 = Iend,)
!  gradv_eb(ie1:ie,1:je)
!-------------------------------------------------------------------------
      if (TanEast) then
      if (EAST .and. have_g_ie) then
         gradv_eb=0.0
         voe_eb0=VVV(ie2:ie,:,:)
         do k=1,ke
            do j=je,J_start+1,-1
               if (AMSUE(IE,J,K) .eq. 1) THEN
               	! dVdy
               gradv_eb(ie1,j) = voe_eb1(ie1,j-1,k) - voe_eb1(ie1,j,k)
               gradv_eb(ie, j) = voe_eb1(ie, j-1,k) - voe_eb1(ie, j,k)
               end if
            end do
! Attentions: the corners excluded
               do j=je2,J_start+2,-1
                  if (AMSUE(IE,J,K) .eq. 1) THEN
                  dVdt=voe_eb1(ie1,j,k)-voe_eb0(ie1,j,k)
                  dVdx=voe_eb0(ie1,j,k)-voe_eb0(ie2,j,k)

                  if ((dVdt*dVdx).lt.0.) dVdt=0.
                  if ((dVdt*(gradv_eb(ie1,j+1)+gradv_eb(ie1,j))).gt.0.)then
                     dVde=gradv_eb(ie1,j+1)
                  else
                     dVde=gradv_eb(ie1,j  )
                  end if
                     cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
                     Cx=dVdt*dVdx
#ifdef RADIATION_2D
                     Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
                     Ce=0.
#endif /*RADIATION_2D*/
                     obc_ve(j,k)  = (cff*voe_eb1(ie,j,k)+               &
     &                               Cx *voe_eb0(ie1,j,k)-              &
     &                               MAX(Ce,0.)*gradv_eb(ie,j+1)-        &
     &                               MIN(Ce,0.)*gradv_eb(ie,j  ))/       &
     &                               (cff+Cx)

#ifdef RADIATION_NUDGING
               if ((dVdt*dVdx) .lt. 0.0) then
                  tau = M3obc_in
               else
                  tau = M3obc_out
               end if
#ifdef IMPLICIT_NUDGING
               if (tau.gt.0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_ve(j,k) = (1.0 - phi)*obc_ve(j,k) +                  &
     &                       phi * (obcf_ve(j,k) + VSET(ie,j))
#else
               obc_ve(j,k) = obc_ve(j,k) +                              &
     &                       tau*((obcf_ve(j,k)+VSET(ie,j))-            &
     &                            voe_eb1(ie,j,k) )
#endif

#endif /*RADIATION_NUDGING*/
               end if
               end do
         end do
         voe_eb1=voe_eb0
      end if
      endif

! Update OBC predict uv to model fields

      CALL UPDATE_ORLANSKI_UV(UUU,VVV) ! in mo_obcs.f90

      UUU0(:,:,:)=UUU(:,:,:)
      VVV0(:,:,:)=VVV(:,:,:)

#endif /*ORLANSKI*/

      END SUBROUTINE ORLANSKI_UV

! -------------------- ORLANSKI TS OBCS ----------------------!
      SUBROUTINE ORLANSKI_TS(TTT0,SSS0,i_step)
#ifdef ORLANSKI
      REAL TTT(IE,JE,KE),TTT0(IE,JE,KE),SSS(IE,JE,KE),SSS0(IE,JE,KE)
      REAL tobc_in, tobc_out
      phi = 0.0
      tau = 0.0
      ! tobc_in = 4.0 * ( 1.0/(1.0*86400.0) ) !==OBC_fac*(1/OBCCOF(DAYS)*86400)
      ! OBCcof: Nudging/Relaxtion time scale. TNUDG_ROMS set to
      ! 2*2.0(days).
      ! OBCfac: Factor between passive(outflow) and active(inflow) open
      !         boundary conditions, If OBCFAC > 1, nudging on inflow is
      !         stronger
      !         than on outflow (which is recommended).
      ! tobc_out = ( 1.0/(1.0*86400.0) )

      tobc_in   = ( 4.0 *1.0/(1.0*86400.0) )
      tobc_out = ( 1.0/(1.0*86400.0) )

      TTT(:,:,:)=TTT0(:,:,:)
      SSS(:,:,:)=SSS0(:,:,:)
! Western edge, all BCs in the ROMS code.
      if (WEST .and. have_g_is) then
       grad_wb=0.0
       tho_wb0 = TTT(1:3,:,:)
       do k=1,ke
         do j=je1,1,-1 ! ATTEN: The mpiom's je is from north to south,oppsite to Roms.
         	 ! dTdy
            grad_wb(1,j) = tho_wb1(1,j,k) - tho_wb1(1,j+1,k)
            grad_wb(2,j) = tho_wb1(2,j,k) - tho_wb1(2,j+1,k)
         enddo
! Attentions: the corners excluded
            do j=je1,2,-1
               dTdt=tho_wb1(2,j,k) - tho_wb0(2,j,k)
               dTdx=tho_wb0(2,j,k) - tho_wb0(3,j,k)   ! upstream scheme

               if ((dTdt*dTdx) .lt. 0.) dTdt=0.
               if ((dTdt*(grad_wb(2,j)+grad_wb(2,j-1))) .gt. 0.) then
                  dTde = grad_wb(2,j)
               else
                  dTde = grad_wb(2,j-1)
               end if
               cff = max(dTdx*dTdx+dTde*dTde,obcerr)
               Cx = dTdt*dTdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dTdt*dTde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
               obc_tw(j,k) = (cff*tho_wb1(1,j,k)+                       &
     &                        Cx *tho_wb0(2,j,k)-                       &
     &                        MAX(Ce,0.)*grad_wb(1,j  )-                &
     &                        MIN(Ce,0.)*grad_wb(1,j-1))/               &
     &                        (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dTdt*dTdx) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_tw(j,k) = obc_tw(j,k) +                              &
     &                       tau*(cli_tw(j,k)-                          &
     &                            tho_wb1(1,j,k) )

#endif /*RADIATION_NUDGING*/
            enddo
       enddo
       tho_wb1=tho_wb0
      end if

!---------------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.(ie1 = Iend, je1
!  =Jstr)
!---------------------------------------------------------------------------
!
      if (EAST .and. have_g_ie) then
         grad_eb=0.0
         tho_eb0=TTT(ie2:ie,:,:)


         do k=1,ke
            do j=je1,1,-1
               grad_eb(ie1,j) = tho_eb1(ie1,j,k) - tho_eb1(ie1,j+1,k)
               grad_eb(ie, j) = tho_eb1(ie, j,k) - tho_eb1(ie ,j+1,k)
            end do
! Attentions: the corners excluded
            do j=je1,2,-1
               dTdt=tho_eb1(ie1,j,k) - tho_eb0(ie1,j,k)
               dTdx=tho_eb0(ie1,j,k) - tho_eb0(ie2,j,k) ! upstream scheme

            if ((dTdt*dTdx) .lt. 0.) dTdt=0.
            if ((dTdt*(grad_eb(ie1,j)+grad_eb(ie1,j-1))) .gt. 0.) then
               dTde = grad_eb(ie1,j)
            else
               dTde = grad_eb(ie1,j-1)
            end if
            cff = max(dTdx*dTdx+dTde*dTde,obcerr)
            Cx = dTdt*dTdx
#ifdef RADIATION_2D
            Ce = MIN(cff,MAX(dTdt*dTde,-cff))
#else
            Ce = 0.
#endif /*RADIATION_2D*/
            obc_te(j,k) = (cff*tho_eb1(ie, j,k)+                        &
     &                     Cx *tho_eb0(ie1,j,k)-                        &
     &                     MAX(Ce,0.)*grad_eb(ie,j  )-                  &
     &                     MIN(Ce,0.)*grad_eb(ie,j-1))/                 &
     &                     (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dTdt*dTdx) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_te(j,k) = obc_te(j,k) +                              &
     &                       tau*(cli_te(j,k)-                          &
     &                            tho_eb1(ie,j,k) )

#endif /*RADIATION_NUDGING*/
            end do
         end do
         tho_eb1=tho_eb0
      endif

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  Later boundary conditions at the southern edge.(1 = I_start = IstrU,
!  ie1 = Iend, je-1 = Jstr)
!----------------------------------------------------------------------------
      if (SOUTH .and. have_g_je) then
         grad_sb=0.0
         tho_sb0=TTT(:,je2:je,:)


         do k=1,ke
            do i=2,ie
               grad_sb(i,je1) = tho_sb1(i,je1 ,k) - tho_sb1(i-1,je1 ,k)
               grad_sb(i,je)  = tho_sb1(i,je,k) - tho_sb1(i-1,je,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
                  dTdt=tho_sb1(i,je1,k) - tho_sb0(i,je1,k)
                  dTde=tho_sb0(i,je1,k) - tho_sb0(i,je2,k) ! upstream scheme

                  if ((dTdt*dTde) .lt. 0.) dTdt=0.
                  if ((dTdt*(grad_sb(i,je1)+grad_sb(i+1,je1))) .gt. 0.)then
                     dTdx = grad_sb(i,je1)
                  else
                     dTdx = grad_sb(i+1,je1)
                  end if
                  cff=MAX(dTdx*dTdx+dTde*dTde,obcerr)
                  Ce=dTdt*dTde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dTdt*dTdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_ts(i,k) = (cff*tho_sb1(i,je, k)+                  &
     &                           Ce *tho_sb0(i,je1,k)-                  &
     &                           MAX(Cx,0.)*grad_sb(i,je)-              &
     &                           MIN(Cx,0.)*grad_sb(i+1,je))/           &
     &                           (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dTdt*dTde) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_ts(i,k) = obc_ts(i,k) +                              &
     &                       tau*(cli_ts(i,k)-                          &
     &                            tho_sb1(i,je,k) )

#endif /*RADIATION_NUDGING*/
               end do
         end do
         tho_sb1=tho_sb0
      endif

!
!-----------------------------------------------------------------------------
!  Later boundary conditions at the northern edge.(1 = I_start = IstrU,
!  2 = Jend)
!-----------------------------------------------------------------------------

      if (NORTH .and. have_g_js) then
         grad_nb=0.0
         tho_nb0=TTT(:,1:3,:)


         do k=1,ke
            do i=2,ie
               grad_nb(i,2) = tho_nb1(i,2,k) - tho_nb1(i-1,2,k)
               grad_nb(i,1) = tho_nb1(i,1,k) - tho_nb1(i-1,1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
                  dTdt=tho_nb1(i,2,k) - tho_nb0(i,2,k)
                  dTde=tho_nb0(i,2,k) - tho_nb0(i,3,k) ! upstream scheme

                  if ((dTdt*dTde) .lt. 0.) dTdt=0.
                  if ((dTdt*(grad_nb(i,2)+grad_nb(i+1,2))) .gt. 0.) then
                     dTdx = grad_nb(i,2)
                  else
                     dTdx = grad_nb(i+1,2)
                  end if
                  cff=MAX(dTdx*dTdx+dTde*dTde,obcerr)
                  Ce=dTdt*dTde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dTdt*dTdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_tn(i,k) = (cff*tho_nb1(i,1,k)+                    &
     &                          Ce *tho_nb0(i,2,k)-                     &
     &                          MAX(Cx,0.)*grad_nb(i  ,1)-              &
     &                          MIN(Cx,0.)*grad_nb(i+1,1))/             &
     &                          (cff+Ce)

#ifdef RADIATION_NUDGING
               if ((dTdt*dTde) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_tn(i,k) = obc_tn(i,k) +                              &
     &                       tau*(cli_tn(i,k)-                          &
     &                            tho_nb1(i,1,k) )

#endif /*RADIATION_NUDGING*/
               end do
         end do
         tho_nb1=tho_nb0
      endif

!---------------------------------------------------------------------------------------------
!========================== Sanlity ============================!
!---------------------------------------------------------------------------------------------

! Western edge, all BCs in the ROMS code.
      if (WEST .and. have_g_is) then
       grad_wb=0.0
       sao_wb0 = SSS(1:3,:,:)


       do k=1,ke
         do j=je1,1,-1 ! ATTEN: The mpiom's je is from north to south,oppsite to Roms.
            grad_wb(1,j) = sao_wb1(1,j,k) - sao_wb1(1,j+1,k)
            grad_wb(2,j) = sao_wb1(2,j,k) - sao_wb1(2,j+1,k)
         enddo
! Attentions: the corners excluded
            do j=je1,2,-1
               dSdt=sao_wb1(2,j,k) - sao_wb0(2,j,k)
               dSdx=sao_wb0(2,j,k) - sao_wb0(3,j,k)

               if ((dSdt*dSdx) .lt. 0.) dSdt=0.
               if ((dSdt*(grad_wb(2,j)+grad_wb(2,j-1))) .gt. 0.) then
                  dSde = grad_wb(2,j)
               else
                  dSde = grad_wb(2,j-1)
               end if
               cff = max(dSdx*dSdx+dSde*dSde,obcerr)
               Cx = dSdt*dSdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dSdt*dSde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
               obc_sw(j,k) = (cff*sao_wb1(1,j,k)+                       &
     &                        Cx *sao_wb0(2,j,k)-                       &
     &                        MAX(Ce,0.)*grad_wb(1,j  )-                &
     &                        MIN(Ce,0.)*grad_wb(1,j-1))/               &
     &                        (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dSdt*dSdx) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_sw(j,k) = obc_sw(j,k) +                              &
     &                       tau*(cli_sw(j,k)-                          &
     &                            sao_wb1(1,j,k) )

#endif /*RADIATION_NUDGING*/
            enddo
       enddo
       sao_wb1=sao_wb0
      end if

!
!---------------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.(ie1 = Iend, je1
!  =Jstr)
!---------------------------------------------------------------------------
!
      if (EAST .and. have_g_ie) then
         grad_eb=0.0
         sao_eb0=SSS(ie2:ie,:,:)


         do k=1,ke
            do j=je1,1,-1
               grad_eb(ie1,j) = sao_eb1(ie1,j,k) - sao_eb1(ie1,j+1,k)
               grad_eb(ie, j) = sao_eb1(ie, j,k) - sao_eb1(ie ,j+1,k)
            end do
! Attentions: the corners excluded
            do j=je1,2,-1
               dSdt=sao_eb1(ie1,j,k) - sao_eb0(ie1,j,k)
               dSdx=sao_eb0(ie1,j,k) - sao_eb0(ie2,j,k)

            if ((dSdt*dSdx) .lt. 0.) dSdt=0.
            if ((dSdt*(grad_eb(ie1,j)+grad_eb(ie1,j-1))) .gt. 0.) then
               dSde = grad_eb(ie1,j)
            else
               dSde = grad_eb(ie1,j-1)
            end if
            cff = max(dSdx*dSdx+dSde*dSde,obcerr)
            Cx = dSdt*dSdx
#ifdef RADIATION_2D
            Ce = MIN(cff,MAX(dSdt*dSde,-cff))
#else
            Ce = 0.
#endif /*RADIATION_2D*/
            obc_se(j,k) = (cff*sao_eb1(ie, j,k)+                        &
     &                     Cx *sao_eb0(ie1,j,k)-                        &
     &                     MAX(Ce,0.)*grad_eb(ie,j  )-                  &
     &                     MIN(Ce,0.)*grad_eb(ie,j-1))/                 &
     &                     (cff+Cx)
#ifdef RADIATION_NUDGING
               if ((dSdt*dSdx) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_se(j,k) = obc_se(j,k) +                              &
     &                       tau*(cli_se(j,k)-                          &
     &                            sao_eb1(ie,j,k) )

#endif /*RADIATION_NUDGING*/
            end do
         end do
         sao_eb1=sao_eb0
      endif

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  Later boundary conditions at the southern edge.(1 = I_start = IstrU,
!  ie1 = Iend, je-1 = Jstr)
!----------------------------------------------------------------------------
      if (SOUTH .and. have_g_je) then
         grad_sb=0.0
         sao_sb0=SSS(:,je2:je,:)


         do k=1,ke
            do i=2,ie
               grad_sb(i,je1)= sao_sb1(i,je1,k) - sao_sb1(i-1,je1,k)
               grad_sb(i,je) = sao_sb1(i,je,k) - sao_sb1(i-1,je,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
                  dSdt=sao_sb1(i,je1,k) - sao_sb0(i,je1,k)
                  dSde=sao_sb0(i,je1,k) - sao_sb0(i,je2,k)

                  if ((dSdt*dSde) .lt. 0.) dSdt=0.
                  if ((dSdt*(grad_sb(i,je1)+grad_sb(i+1,je1))) .gt.0.)then
                     dSdx = grad_sb(i,je1)
                  else
                     dSdx = grad_sb(i+1,je1)
                  end if
                  cff=MAX(dSdx*dSdx+dSde*dSde,obcerr)
                  Ce=dSdt*dSde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dSdt*dSdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_ss(i,k) = (cff*sao_sb1(i,je, k)+                  &
     &                           Ce *sao_sb0(i,je1,k)-                  &
     &                           MAX(Cx,0.)*grad_sb(i,je)-              &
     &                           MIN(Cx,0.)*grad_sb(i+1,je))/           &
     &                           (cff+Ce)
#ifdef RADIATION_NUDGING
               if ((dSdt*dSde) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_ss(i,k) = obc_ss(i,k) +                              &
     &                       tau*(cli_ss(i,k)-                          &
     &                            sao_sb1(i,je,k) )

#endif /*RADIATION_NUDGING*/
               end do
         end do
         sao_sb1=sao_sb0
      endif

!-----------------------------------------------------------------------------
!  Later boundary conditions at the northern edge.(1 = I_start = IstrU,
!  2 = Jend)
!-----------------------------------------------------------------------------

      if (NORTH .and. have_g_js) then
         grad_nb=0.0
         sao_nb0=SSS(:,1:3,:)


         do k=1,ke
            do i=2,ie
               grad_nb(i,2) = sao_nb1(i,2,k) - sao_nb1(i-1,2,k)
               grad_nb(i,1) = sao_nb1(i,1,k) - sao_nb1(i-1,1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
                  dSdt=sao_nb1(i,2,k) - sao_nb0(i,2,k)
                  dSde=sao_nb0(i,2,k) - sao_nb0(i,3,k)

                  if ((dSdt*dSde) .lt. 0.) dSdt=0.
                  if ((dSdt*(grad_nb(i,2)+grad_nb(i+1,2))) .gt. 0.) then
                     dSdx = grad_nb(i,2)
                  else
                     dSdx = grad_nb(i+1,2)
                  end if
                  cff=MAX(dSdx*dSdx+dSde*dSde,obcerr)
                  Ce=dSdt*dSde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dSdt*dSdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_sn(i,k) = (cff*sao_nb1(i,1,k)+                    &
     &                          Ce *sao_nb0(i,2,k)-                     &
     &                          MAX(Cx,0.)*grad_nb(i,1)-                &
     &                          MIN(Cx,0.)*grad_nb(i+1,1))/             &
     &                          (cff+Ce)
#ifdef RADIATION_NUDGING
               if ((dSdt*dSde) .lt. 0.0) then
                  tau = tobc_in
               else
                  tau = tobc_out
               end if
               tau = tau*dt

               obc_sn(i,k) = obc_sn(i,k) +                              &
     &                       tau*(cli_sn(i,k)-                          &
     &                            sao_nb1(i,1,k) )

#endif /*RADIATION_NUDGING*/
               end do
         end do
         sao_nb1=sao_nb0
      endif

      CALL UPDATE_ORLANSKI_TS(TTT,SSS) ! in mo_obcs.f90

      TTT0(:,:,:)=TTT(:,:,:)
      SSS0(:,:,:)=SSS(:,:,:)

#endif /*ORLANSKI*/

      END SUBROUTINE ORLANSKI_TS

      SUBROUTINE  OBCS_TST_UVO(I_step)
! Tidal and Subtidal Open Boundary Routine of barotropic component
! H.H 23.5.20 upgraded this code for nesting simulation
#ifdef OBCS_TST_NEST
      REAL M2FC,M2EtC
      REAL M2obc_in, M2obc_out, M2nudgcof
      INTEGER inflow,I_step, obcfac

      obcfac = 4.0
      M2nudgcof = (2.0/24.0)*86400
      M2obc_out = M2nudgcof
      M2obc_in =obcfac * M2obc_out

      M2FC  = 1.0
      M2EtC = 1.0
      tau = 0.0
! -------------------------------- U_bt ----------------------------------------------
! Western Edge
      if (WEST .and. have_g_is .and. TSTWest) then
      	gradu_wb = 0.0
      	! U_Rb = Ubar - U_Tb
      	USOR_wb0(1:3,:) = USO(1:3,:) - USOT(1:3,:)  ! kout
         if (I_STEP .eq. 1 ) then
           ! USOR_wb1 = USOR_wb0
! Update the OB values (only index=1)
!            do j=je,1,-1 ! Corners included
!              if (WETO(1,J,1) .eq. 1) THEN
!                cff = 1.0 / (0.5*(DEPTO(1,j)+ZOOLD(1,j) + DEPTO(2,j)+ ZOOLD(2,j) ) )
!                Cx  = SQRT(G*cff)*M2FC
!                elevation_tiaozheng =-1.0*Cx*( 0.5* (ZOOLD(1,j) + ZOOLD(2,j) )      &
!     &                           -obc_zTNw(j)-obc_zsubNw(j)*M2EtC)
!              !  USOR_wb0(1,j) = USOR_wb0(1,j) - elevation_tiaozheng
!                USOR_wb1(1,j) = USOR_wb1(1,j) - elevation_tiaozheng ! know
!              end if
!            enddo
          endif
!  U_Tb = USOT +/- elevation_tiaozheng
!  U_Rb = Ubar - U_Tb
! Residual components of the OB values only have been left
! namely, USOR_wb0(1,1:je) & USOR_wb1(1,1:je) are U_Rb on OB
!Current step: USOR_wb0(1:3,1:je), obc_ubarw,   zeta_subT,  zetaT
!Former step: USOR_wb1(1:3,1:je), obc_ubarNw, obc_zsubNw, obc_zTNw
!Current step: USOR_eb0(ie-3:ie1,1:je), obc_ubarw,   zeta_subT,  zetaT
!Former step: USOR_eb1(ie-3:ie1,1:je), obc_ubarNw, obc_zsubNe, obc_zTNe

! So 2-D Radiation Scheme for U_Rgb can start
! U_Rgb = U_Rb - U_Rlb = Ubar - U_Tb - U_Rlb
! dUdy where means d(U_Rgb) / dy
         do j=je1,1,-1    ! ROMS: DO j=Jstr-1,Jend+1
            if (AMSUO(0,J,1) .eq. 1) THEN
            gradu_wb(0,j) = USOR_wb1(1,j)   - obc_ubarNw(j) -            &
     &                     USOR_wb1(1,j+1) + obc_ubarNw(j+1)  ! OB
            gradu_wb(1,j) = USOR_wb1(2,j)   - obc_ubarNw(j) -            &
     &                     USOR_wb1(2,j+1) + obc_ubarNw(j+1)  ! OB+1
            end if
         end do
! Attenetion:  the corners excluded
        do j=je1,2,-1 ! ROMS: DO j=Jstr,Jend+1 ubar(Istr,j,kout)
! dUdt & dUdx
! d(U_Rgb) / dt & d(U_Rgb) / dx
          ! WRITE(IO_STDOUT,*) '  West TST-BT Main begin ', i
           if (AMSUO(0,J,1) .eq. 1) THEN
           dUdt=USOR_wb1(2,j) - obc_ubarNw(j) -                         &
     &          (USOR_wb0(2,j) - obc_ubarw(j))
           dUdx=USOR_wb0(2,j) - obc_ubarw(j)  -                         &
     &          (USOR_wb0(3,j) - obc_ubarw(j))   ! upstream shceme

           if ((dUdt*dUdx) .lt. 0.0) then
              inflow = 1
              tau = M2obc_in
           else
              tau = M2obc_out
              inflow = 0
           end if

#ifdef IMPLICIT_NUDGING
           IF (tau .gt. 0.0) tau=1.0/tau
#else
           tau = tau*dt
#endif
           if ((dUdt*dUdx) .lt. 0.) dUdt=0.
           if ((dUdt*(gradu_wb(1,j)+gradu_wb(1,j-1))) .gt. 0.) then ! upstream shceme
              dUde = gradu_wb(1,j)  ! OB+1
           else
              dUde = gradu_wb(1,j-1)
           end if
           cff = max(dUdx*dUdx+dUde*dUde,obcerr)
           Cx = dUdt*dUdx
#ifdef RADIATION_2D
           Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
           Ce = 0.
#endif /*RADIATION_2D*/

           USO(1,j) = (cff* (USOR_wb1(1,j) - obc_ubarNw(j)) +           &
     &                  Cx* (USOR_wb0(2,j) - obc_ubarw(j) ) -           &
     &                        MAX(Ce,0.)*gradu_wb(0,j  )-                &
     &                        MIN(Ce,0.)*gradu_wb(0,j-1))/               &
     &                        (cff+Cx)

           cff = 1.0 / (0.5*(DEPTO(1,j)+ZOOLD(1,j) + DEPTO(2,j)+ ZOOLD(2,j) ) )
           Cx  = SQRT(G*cff)*M2FC
           elevation_tiaozheng =-1.0*Cx*( 0.5*(ZOOLD(1,j) - zetaT(1,j)       &
     &                          +ZOOLD(2,j) - zetaT(2,j))                                        &
     &                          - zeta_subT(1,j)*M2EtC )

           if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
           phi = dt / (tau + dt)
           USO(1,j) = USOR_wb1(1,j) -                                   &
     &                phi*(USOR_wb1(1,j) - obc_ubarNw(j)) +             &
     &                USOT(1,j) + elevation_tiaozheng
#else
           USO(1,j) = USOR_wb1(1,j) -                                   &
     &                tau*(USOR_wb1(1,j) - obc_ubarNw(j)) +             &
     &                USOT(1,j) + elevation_tiaozheng
#endif
           else
           USO(1,j) = USO(1,j)  + obc_ubarw(j) +                        &
     &                USOT(1,j) + elevation_tiaozheng
           endif
           USOR_wb0(1,j) = USO(1,j) - USOT(1,j) - elevation_tiaozheng ! OB values updated
           obc_ubarNw(j) = obc_ubarw(j)  ! subtidal current on OB of the former step
           obc_zsubNw(j)    = zeta_subT(1,j) ! subtidal elevation on OB of the former step
           obc_zTNw(j)   =0.5*(zetaT(1,j) + zetaT(2,j)) ! Tidal elevation on OB
! Masking
           USO(1,j) = USO(1,j) * amsuo(1,j,1)
           USOR_wb0(1,j) = USOR_wb0(1,j) * amsuo(1,j,1)
           USOR_wb1(:,j) = USOR_wb0(:,j)
           end if
         end do
! Southwest Corner
         obc_ubarNw(je) = obc_ubarw(je)
         obc_zTNw(je)   = 0.5*(zetaT(1,je) + zetaT(2,je) )
         obc_zsubNw(je)    = zeta_subT(1,je)  ! obc_zw(je)
! Northwest Corner
         obc_ubarNw(1)  = obc_ubarw(1)
         obc_zTNw(1)    = 0.5*(zetaT(1,1) + zetaT(2,1) )
         obc_zsubNw(1)     = zeta_subT(1,1) ! obc_zw(1)
         ! WRITE(IO_STDOUT,*) ' West TST-BT ', I_STEP
         ! WRITE(IO_STDOUT,*) ' West TST:  USOR_wb0 (U_Rb) ', USOR_wb0(1,:)
         ! WRITE(IO_STDOUT,*) ' West TST:  USOT (U_Tlb) ', USOT(1,:)
         ! WRITE(IO_STDOUT,*) ' West TST:  obc_ubarw (U_Rlb) ', obc_ubarw(:)
      end if

! Eastern Edge
      if (East .and. have_g_ie .and. TSTEast ) then
         gradu_eb = 0.0
         USOR_eb0(ie-3:ie1,:)= USO(ie-3:ie1,:) - USOT(ie-3:ie1,:) ! kout
         if (I_STEP .eq. 1 ) then
           ! USOR_eb1 = USOR_eb0
!           do j=je,1,-1  ! ROMS: DO j=Jstr-1,Jend+1
!            if (WETO(IE,J,1) .eq. 1) THEN
!            cff = 1.0 / (0.5*(DEPTO(ie,j) + ZOOLD(ie,j)+DEPTO(ie1,j) + ZOOLD(ie1,j)))
!            Cx  = SQRT(G*cff)*M2FC
!            elevation_tiaozheng =Cx*(0.5* (ZOOLD(ie,j) + ZOOLD(ie1,j))             &
!     &                           -obc_zTNe(j)-obc_zsubNe(j)*M2EtC)
           !  USOR_eb0(ie1,j) = USOR_eb0(ie1,j) - elevation_tiaozheng
!            USOR_eb1(ie1,j) = USOR_eb1(ie1,j) - elevation_tiaozheng ! know
!            endif
!           enddo
         endif
! dUdy where means d(U_Rgb = U_Rb - U_Rlb) / dy
         do j=je1,1,-1 ! ROMS: DO j=Jstr,Jend+1  ubar(Iend+1,j,kout)
            if (AMSUO(IE,J,1) .eq. 1) THEN
            gradu_eb(ie1,j) = USOR_eb1(ie2,j)   - obc_ubarNe(j) -        &
     &                     USOR_eb1(ie2,j+1) + obc_ubarNe(j+1) ! OB+1
            gradu_eb(ie ,j) = USOR_eb1(ie1,j)   - obc_ubarNe(j) -         &
     &                     USOR_eb1(ie1,j+1) + obc_ubarNe(j+1) ! OB
            end if
         end do
! Attenetion:  the corners excluded
         do j=je1,2,-1
           if (AMSUO(IE,J,1) .eq. 1) THEN
           dUdt=USOR_eb1(ie2,j) - obc_ubarNe(j) -                       &
     &          USOR_eb0(ie2,j) + obc_ubare(j)
           dUdx=USOR_eb0(ie2,j) - obc_ubare(j)  -                       &
     &          USOR_eb0(ie-3,j) + obc_ubare(j)   ! upstream shceme

           if ((dUdt*dUdx) .lt. 0.0) then
              inflow = 1
              tau = M2obc_in
           else
              tau = M2obc_out
              inflow = 0
           end if
#ifdef IMPLICIT_NUDGING
           IF (tau .gt. 0.0) tau=1.0/tau
#else
           tau = tau*dt
#endif
         if ((dUdt*dUdx) .lt. 0.) dUdt=0.
           if ((dUdt*(gradu_eb(ie1,j)+gradu_eb(ie1,j-1))) .gt. 0.) then ! upstream shceme
              dUde = gradu_eb(ie1,j)  ! OB+1
           else
              dUde = gradu_eb(ie1,j-1)
           end if
           cff = max(dUdx*dUdx+dUde*dUde,obcerr)
           Cx = dUdt*dUdx
#ifdef RADIATION_2D
           Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
           Ce = 0.
#endif /*RADIATION_2D*/

         USO(ie1,j) = (cff* (USOR_eb1(ie1,j)  - obc_ubarNe(j)) +          &
     &                 Cx* (USOR_eb0(ie2,j) - obc_ubare(j) ) -          &
     &                        MAX(Ce,0.)*gradu_eb(ie,j  )-               &
     &                        MIN(Ce,0.)*gradu_eb(ie,j-1))/              &
     &                        (cff+Cx)

           cff = 1.0 /  (0.5*(DEPTO(ie,j)+ZOOLD(ie,j) + DEPTO(ie1,j)+ ZOOLD(ie1,j) ) )
           Cx  = SQRT(G*cff)*M2FC
           elevation_tiaozheng =Cx*( 0.5*(ZOOLD(ie,j) - zetaT(ie,j)          &
     &                          +ZOOLD(ie1,j) - zetaT(ie1,j))                               &
     &                          -zeta_subT(ie,j)*M2EtC )

           if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
           phi = dt / (tau + dt)
           USO(ie1,j) = USOR_eb1(ie1,j) -                                 &
     &                 phi*(USOR_eb1(ie1,j) - obc_ubarNe(j)) +           &
     &                 USOT(ie1,j) + elevation_tiaozheng
#else
           USO(ie1,j) = USOR_eb1(ie1,j) -                                 &
     &                tau*(USOR_eb1(ie1,j) - obc_ubarNe(j)) +            &
     &                USOT(ie1,j) + elevation_tiaozheng
#endif
           else
           USO(ie1,j) = USO(ie1,j)  + obc_ubare(j) +                      &
     &                 USOT(ie1,j) + elevation_tiaozheng
           endif
           USOR_eb0(ie1,j) = USO(ie1,j) - USOT(ie1,j) - elevation_tiaozheng ! OB values updated
           obc_ubarNe(j)  = obc_ubare(j)  ! SubTidal Current on OB
           obc_zsubNe(j)     = zeta_subT(ie,j) ! SubTidal elevation on OB
           obc_zTNe(j)    = 0.5*(zetaT(ie,j)+zetaT(ie1,j)) ! Tidal elevation on OB
! MASKING
           USO(ie1,j) = USO(ie1,j) * amsuo(ie1,j,1)
           USOR_eb0(ie1,j) = USOR_eb0(ie1,j) * amsuo(ie1,j,1)
           USOR_eb1(:,j)  = USOR_eb0(:,j)
           end if
         end do
! SouthEast Corner
           obc_ubarNe(je) = obc_ubare(je)
           obc_zTNe(je)   = 0.5*(zetaT(ie,je) + zetaT(ie1,je) )
           obc_zsubNe(je)    = zeta_subT(ie,je)
! NorthEast Corner
           obc_ubarNe(1)  = obc_ubare(1)
           obc_zTNe(1)    =  0.5*(zetaT(ie,1) + zetaT(ie1,1) )
           obc_zsubNe(1)     = zeta_subT(ie,1)
          ! WRITE(IO_STDOUT,*) ' East TST-BT ', I_STEP
      endif


! Southern Edge
! ATTEN: This is also Chapman BC.
        if (ChpmSouth .and.  TSTSouth) then
        if (SOUTH .and. have_g_je )  then
            if (WEST .and. have_g_is) then
        	     i_begin = 2 ! SouthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1
           	   i_end = ie2  ! SouthEast Corner Excluded
           	else
               i_begin = 1
           	   i_end = ie1
           endif

           do i=i_begin,i_end
              if (AMSUO(I,JE,1) .eq. 1) THEN
                 cff  = dpyo(i,je1)  ! DPYO = DT / DLYU
              ! ATTEN: JE1 is for NORTH - SOUTH.
              ! C = - SQRT(G*D)
      !          cff1 = SQRT(G*(DEUTO(i,je1) +                           &
      ! &               0.5 * (ZOOLD(i,je1) + ZOOLD(i+1,je1))))
                 cff1 = SQRT(G*0.5*(DEPTO(i,je1) +ZOOLD(i,je1)     &
     &                     +DEPTO(i+1,je1)+ ZOOLD(i+1,je1)))
                 Ce   = cff * cff1
                 cff2 = 1.0 / (1.0 + Ce)
                 USO(i,je) = cff2 * (USOOLD(i,je) +                          &
     &                         Ce*USO(i,je1))
                 USO(i,je) = USO(i,je) * AMSUO(i,je,1)

              end if
          end do
         endif
         endif

! Northern Edge
! ATTEN: This is also Chapman BC.
       if (ChpmNorth .and.  TSTNorth ) then
       if (NORTH .and. have_g_js) then
            if (WEST .and. have_g_is) then
        	     i_begin = 2  ! NorthWest Corner Excluded
        	     i_end = ie
           else if (East .and. have_g_ie) then
           	   i_begin = 1
           	   i_end = ie2   ! NorthEast Corner Excluded
           	else
               i_begin = 1
           	   i_end = ie1
           endif
           do i=i_begin,i_end
              if (AMSUO(I,1,1) .eq. 1) THEN
                cff = dpyo(i,2) ! DPYO = DT / DLYU
                ! ATTEN: 2 is for NORTH - SOUTH.
                ! C = SQRT(G*D)
!               cff1 = SQRT(G*(DEUTO(i,2) +                               &
!    &               0.5 * (ZOOLD(i,2) + ZOOLD(i+1,2))))
                cff1 = SQRT(G*0.5*(DEPTO(i,2) +ZOOLD(i,2)     &
     &                     +DEPTO(i+1,2)+ ZOOLD(i+1,2)))
                Ce   = cff * cff1
                cff2 = 1.0 / (1.0 + Ce)
                USO(i,1) = cff2 * (USOOLD(i,1) +                           &
     &                         Ce*USO(i,2))
                USO(i,1) = USO(i,1) * AMSUO(i,1,1)
            end if
         end do
        endif
       endif

! -------------------------------- V_bt ----------------------------------------------
! Southern Edge
      if (SOUTH .and. have_g_je .and.  TSTSouth) then
      	 gradv_sb(:,:) = 0.0
         VSER_sb0(:,je-3:je1) = VSE(:,je-3:je1) - VSET(:,je-3:je1)
         if (I_STEP .eq. 1 ) then
          ! VSER_sb1 = VSER_sb0
!         do i=1,ie
!               if (WETO(I,JE,1) .eq. 1) THEN
!               cff = 1.0 / (0.5*(DEPTO(i,je) + ZOOLD(i,je)+DEPTO(i,je1) + ZOOLD(i,je1)))
!               Ce  = SQRT(G*cff)*M2FC
!               elevation_tiaozheng =-1.0*Ce*(0.5* (ZOOLD(i,je) + ZOOLD(i,je1))       &
!     &                              -obc_zTNs(i)-obc_zsubNs(i)*M2EtC)
!             !  VSER_sb0(i,je) = VSER_sb0(i,je) - elevation_tiaozheng
!               VSER_sb1(i,je) = VSER_sb1(i,je) - elevation_tiaozheng
!               end if
!           end do
         endif
! dVdx
           do i=2,ie
               if (AMSUE(I,JE,1) .eq. 1) THEN
               gradv_sb(i,je)  = VSER_sb1(i,je1)   - obc_vbarNs(i) -         &
     &                       VSER_sb1(i-1,je1) + obc_vbarNs(i-1) ! OB
               gradv_sb(i,je1) = VSER_sb1(i,je2)  - obc_vbarNs(i) -         &
     &                       VSER_sb1(i-1,je2) + obc_vbarNs(i-1) ! OB+1
               endif
           enddo
! Attentions: the corners excluded
           do i=2,ie1
            !  WRITE(IO_STDOUT,*) ' South TST-BT Main begin ', i
               if (AMSUE(I,JE,1) .eq. 1) THEN
                  dVdt=VSER_sb1(i,je2) - obc_vbarNs(i) -                       &
     &                 VSER_sb0(i,je2) + obc_vbars(i)
                  dVde=VSER_sb0(i,je2) - obc_vbars(i)  -                       &
     &                 VSER_sb0(i,je-3) + obc_vbars(i)   ! upstream shceme

             if ((dVdt*dVde) .lt. 0.0) then
                inflow = 1
                tau = M2obc_in
             else
                tau = M2obc_out
                inflow = 0
             end if
#ifdef IMPLICIT_NUDGING
              IF (tau .gt. 0.0) tau=1.0/tau
#else
              tau = tau*dt
#endif
             if ((dVdt*dVde) .lt. 0.) dVdt=0.
             if ((dVdt*(gradv_sb(i,je1)+gradv_sb(i+1,je1))) .gt. 0.) then ! upstream shceme
                dVdx = gradv_sb(i,je1)  ! OB+1
             else
                dVdx = gradv_sb(i+1,je1)  ! OB+1
             end if
             cff = max(dVdx*dVdx+dVde*dVde,obcerr)
             Ce = dVdt*dVde
#ifdef RADIATION_2D
             Cx = MIN(cff,MAX(dVdt*dVdx,-cff))
#else
             Cx = 0.
#endif /*RADIATION_2D*/

             VSE(i,je1) = (cff* (VSER_sb1(i,je1) - obc_vbarNs(i)) +         &
     &                  Ce* (VSER_sb0(i,je2) - obc_vbars(i) ) -         &
     &                        MAX(Cx,0.)*gradv_sb(i,  je )-              &
     &                        MIN(Cx,0.)*gradv_sb(i+1,je))/              &
     &                        (cff+Ce)

             cff = 1.0 /  (0.5*(DEPTO(i,je) + ZOOLD(i,je) +DEPTO(i,je1) + ZOOLD(i,je1)))
             Ce  = SQRT(G*cff)*M2FC
             elevation_tiaozheng =-1.0*Ce*( 0.5*(ZOOLD(i,je) - zetaT(i,je)     &
     &                         +ZOOLD(i,je1) - zetaT(i,je1) )     &
     &                          -zeta_subT(i,je)*M2EtC )

             if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
             phi = dt / (tau + dt)
             VSE(i,je1) = VSER_sb1(i,je1) -                                 &
     &                 phi*(VSER_sb1(i,je1) - obc_vbarNs(i)) +           &
     &                 VSET(i,je1) + elevation_tiaozheng
#else
             VSE(i,je1) = VSER_sb1(i,je1) -                                 &
     &                 tau*(VSER_sb1(i,je1) - obc_vbarNs(i)) +           &
     &                 VSET(i,je1) + elevation_tiaozheng
#endif
             else
             VSE(i,je1) = VSE(i,je1)  + obc_vbars(i) +                      &
     &                 VSET(i,je1) + elevation_tiaozheng
             endif
           VSER_sb0(i,je1) = VSE(i,je1) - VSET(i,je1) - elevation_tiaozheng
           obc_vbarNs(i) = obc_vbars(i)
           obc_zsubNs(i)    = zeta_subT(i,je)
           obc_zTNs(i)   =  0.5*(zetaT(i,je) + zetaT(i,je1) ) ! Tidal elevation on OB
! MASKING
           VSE(i,je1) = VSE(i,je1) * amsue(i,je1,1)
           VSER_sb0(i,je1) = VSER_sb0(i,je1) * amsue(i,je1,1)
           VSER_sb1(i,:)  = VSER_sb0(i,:)
            endif
           enddo
! SouthWest Corner
           obc_vbarNs(1)  = obc_vbars(1)
           obc_zTNs(1)    = zetaT(1,je)
           obc_zsubNs(1)     = zeta_subT(1,je)
! SouthEast Corner
           obc_vbarNs(ie) = obc_vbars(ie)
           obc_zTNs(ie)   = zetaT(ie,je)
           obc_zsubNs(ie)    = zeta_subT(ie,je)
           ! WRITE(IO_STDOUT,*) ' South TST-BT ', I_STEP
          ENDIF

! Northern Edge
          if (NORTH .and. have_g_js .and.  TSTNorth) then
          	gradv_nb(:,:) = 0.0
            VSER_nb0(:,1:3) = VSE(:,1:3) - VSET(:,1:3) ! kout
           if (I_STEP .eq. 1 ) then
            !  VSER_nb1 = VSER_nb0
!            do i=1,ie
!               if (WETO(I,1,1) .eq. 1) THEN
!               cff = 1.0 / (0.5*(DEPTO(i,1) + ZOOLD(i,1) + DEPTO(i,2) + ZOOLD(i,2)  ) )
!               Ce  = SQRT(G*cff)*M2FC
!               elevation_tiaozheng = 1.0*Ce*(0.5*(ZOOLD(i,1)  + ZOOLD(i,2) )              &
!     &                              -obc_zTNn(i)-obc_zsubNn(i)*M2EtC)
!              ! VSER_nb0(i,1) = VSER_nb0(i,1) - elevation_tiaozheng
!               VSER_nb1(i,1) = VSER_nb1(i,1) - elevation_tiaozheng ! know
!               endif
!             enddo
           endif
! dVdx
             do i=2,ie
                if (AMSUE(I,0,1) .eq. 1) THEN
                 gradv_nb(i,1)  = VSER_nb1(i,2)   - obc_vbarNn(i) -           &
     &                      VSER_nb1(i-1,2) + obc_vbarNn(i-1)  ! OB+1
                 gradv_nb(i,0) = VSER_nb1(i,1)  - obc_vbarNn(i) -             &
     &                     VSER_nb1(i-1,1) + obc_vbarNn(i-1) ! OB
                end if
              enddo
! Attentions: the corners excluded
              do i=2,ie1
                !  WRITE(IO_STDOUT,*) ' North TST-BT Main begin ', i
                if (AMSUE(I,0,1) .eq. 1) THEN
                dVdt=VSER_nb1(i,2) - obc_vbarNn(i) -                         &
     &             VSER_nb0(i,2) + obc_vbarn(i)
                dVde=VSER_nb0(i,2) - obc_vbarn(i)  -                         &
     &            VSER_nb0(i,3) + obc_vbarn(i)  ! upstream shceme

                if ((dVdt*dVde) .lt. 0.0) then
                   inflow = 1
                   tau = M2obc_in
                else
                    tau = M2obc_out
                    inflow = 0
                end if
#ifdef IMPLICIT_NUDGING
                IF (tau .gt. 0.0) tau=1.0/tau
#else
                tau = tau*dt
#endif
                if ((dVdt*dVde) .lt. 0.) dVdt=0.
                if ((dVdt*(gradv_nb(i,1)+gradv_nb(i+1,1))) .gt. 0.) then  ! upstream shceme
                    dVdx = gradv_nb(i,1)  ! OB+1
                else
                    dVdx = gradv_nb(i+1,1)  ! OB+1
                end if
                cff = max(dVdx*dVdx+dVde*dVde,obcerr)
                Ce = dVdt*dVde
#ifdef RADIATION_2D
                 Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#else
                 Cx=0.0
#endif /*RADIATION_2D*/
                VSE(i,1) = (cff* (VSER_nb1(i,1) - obc_vbarNn(i)) +           &
     &                     Ce* (VSER_nb0(i,2) - obc_vbarn(i) ) -           &
     &                          MAX(Cx,0.)*gradv_nb(i,  0)    -            &
     &                          MIN(Cx,0.)*gradv_nb(i+1,0))/               &
     &                          (cff+Ce)

                cff = 1.0 /(0.5* (DEPTO(i,1) + ZOOLD(i,1) +DEPTO(i,2) + ZOOLD(i,2) ) )
                Ce  = SQRT(G*cff)*M2FC
                elevation_tiaozheng = Ce*( 0.5*(ZOOLD(i,1) - zetaT(i,1)           &
     &                        +ZOOLD(i,2) - zetaT(i,2)  )   &
     &                          -zeta_subT(i,1)*M2EtC )
                  if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                  phi = dt / (tau + dt)
                  VSE(i,1) = VSER_nb1(i,1) -                                   &
     &                phi*(VSER_nb1(i,1) - obc_vbarNn(i)) +             &
     &                VSET(i,1) + elevation_tiaozheng
                   ! WRITE(IO_STDOUT,*) ' North TST-BT Main nudging ', i
#else
                  VSE(i,1) = VSER_nb1(i,1) -                                   &
     &                tau*(VSER_nb1(i,1) - obc_vbarNn(i)) +             &
     &                VSET(i,0) + elevation_tiaozheng
#endif
                  else
                  VSE(i,1) = VSE(i,1)  + obc_vbarn(i) +                        &
     &                VSET(i,1) + elevation_tiaozheng
                  endif
                 !  WRITE(IO_STDOUT,*) ' North TST-BT Main loop end ', i
                  VSER_nb0(i,1) = VSE(i,1) - VSET(i,1) - elevation_tiaozheng
                  obc_vbarNn(i) = obc_vbarn(i)
                  obc_zsubNn(i)    = zeta_subT(i,1)
                  obc_zTNn(i)   = 0.5*(zetaT(i,1) + zetaT(i,2) ) ! Tidal elevation on OB
! MASKING
                  VSE(i,1) = VSE(i,1) * amsue(i,1,1)
                  VSER_nb0(i,1) = VSER_nb0(i,1) * amsue(i,1,1)
                  VSER_nb1(i,:) = VSER_nb0(i,:)
           end if
         end do
! NorthWest Corner
           obc_vbarNn(1)  = obc_vbarn(1)
           obc_zTNn(1)    = 0.5*(zetaT(1,1) + zetaT(1,2) )
           obc_zsubNn(1)     = zeta_subT(1,1)
! NorthEastCorner
           obc_vbarNn(ie) = obc_vbarn(ie)
           obc_zTNn(ie)   =  0.5*(zetaT(ie,1) + zetaT(ie,2) )
           obc_zsubNn(ie)    =  zeta_subT(1,je)
          ! WRITE(IO_STDOUT,*) ' North TST-BT ', I_STEP
          endif

! Western Edge
      if (ChpmWest .and. TSTWest) then
      if (WEST .and. have_g_is) then
      	if (North .and. have_g_js) then
      		  j_begin = 2  ! NorthWest Corner Excluded
      		  j_end = je
      	else if (South .and. have_g_je) then
      		  j_begin = 1
      		  j_end = je2  ! SouthWest Corner Excluded
        else
      		  j_begin = 1
      		  j_end = je1
      	endif

         do j=j_end,j_begin,-1
           if (AMSUE(1,J,1) .eq. 1) THEN
           cff  = DTDXUE(2,j) ! DTXDUE = DT / DLXV
           ! ATTEN: 2 is for NORTH - SOUTH.
           ! C = - SQRT(G*D)
 !          cff1 = SQRT(G*(DEUTE(2,j) +                               &
 !    &               0.5 * (ZOOLD(2,j ) + ZOOLD(2,j+1))))
           cff1 = SQRT(G*0.5*(DEPTO(2,j) +ZOOLD(2,j)     &
     &                     +DEPTO(2,j+1)+ ZOOLD(2,j+1)))

           Cx   = cff * cff1
           cff2 = 1.0 / (1.0 + Cx)
           VSE(1,j) = cff2 * (VSEOLD(1,j) +                             &
     &                         Cx*VSE(2,j))
           VSE(1,j) = VSE(1,j) * AMSUE(1,j,1)
           end if
        end do
      end if
      endif

! Eastern Edge
      if (ChpmEast .and. TSTEast) then
      if (EAST .and. have_g_ie) then
      	if (North .and. have_g_js) then
      		  j_begin = 2   ! NorthEast Corner Excluded
      		  j_end = je
      	else if (South .and. have_g_je) then
      		  j_begin = 1
      		  j_end = je2    ! SouthEast Corner Excluded
        else
      		  j_begin = 1
      		  j_end = je1
      	endif

         do j=j_end,j_begin,-1
           if (AMSUE(IE,J,1) .eq. 1) THEN
           cff  = DTDXUE(ie1,j) ! DTXDUE = DT / DLXV
           ! ATTEN: JE1 is for NORTH - SOUTH.
           ! C =  SQRT(G*D)
!          cff1 = SQRT(G*(DEUTE(ie1,j) +                             &
!     &               0.5 * (ZOOLD(ie1,j ) + ZOOLD(ie1,j+1))))
           cff1 = SQRT(G*0.5*(DEPTO(ie1,j) +ZOOLD(ie1,j)     &
     &                     +DEPTO(ie1,j+1)+ ZOOLD(ie1,j+1)))

           Cx   = cff * cff1
           cff2 = 1.0 / (1.0 + Cx)
           VSE(ie,j) = cff2 * (VSEOLD(ie,j) +                           &
     &                         Cx*VSE(ie1,j))
           VSE(ie,j) = VSE(ie,j) * AMSUE(ie,j,1)
           end if
         end do
      end if
      endif

      if (Cornersmean) then
! ------------------ BOUNDARY CORNERS BEGIN --------------------------------
! NorthWest Corner
      if (have_g_is .and. have_g_js) then
         if (WEST .and. NORTH) then
         	  if (WETO(1,1,1) .eq. 1) then
              USO(I_start+1,1) = 0.5*(USO(I_start+2,1) +USO(I_start+1,2))*AMSUO(I_start+1,1,1)
              VSE(1,1) = 0.5*(VSE(1,2) +VSE(2,1) )*AMSUE(1,1,1)
           endif
        endif
     endif
! SouthWest Corner
      if (have_g_is .and. have_g_je) then
      	 if (WEST .and. SOUTH) then
      	 	 if (WETO(1,je,1) .eq. 1) then
             USO(I_start+1,je) = 0.5*(USO(I_start+2,je)+USO(I_start+1,je1) )*AMSUO(I_start+1,je,1)
             VSE(1,je1) = 0.5*(VSE(2,je1) +VSE(1,je2) )*AMSUE(1,je1,1)
           endif
        endif
     endif
! NorthEast Corner
      if (have_g_ie .and. have_g_js) then
         if (EAST .and. NORTH) then
         	  if (WETO(ie,1,1) .eq. 1) then
               USO(ie1,1) = 0.5*(USO(ie2,1) +USO(ie1,2) )*AMSUO(ie1,1,1)
               VSE(ie,1) = 0.5*(VSE(ie1,1) +VSE(ie,2) )*AMSUE(ie,1,1)
           endif
        endif
     endif
! SouthEast Corner
      if (have_g_ie .and. have_g_je) then
         if (EAST .and. SOUTH) then
         	  if (WETO(ie,je,1) .eq. 1) then
              USO(ie1,je) = 0.5*(USO(ie2,je) +USO(ie1,je1) )*AMSUO(ie1,je,1)
              VSE(ie,je1) = 0.5*(VSE(ie,je2) +VSE(ie1,je1) )*AMSUE(ie,je1,1)
          endif
        endif
     endif
! ------------------ BOUNDARY CORNERS END --------------------------------
     endif

      CALL bounds_exch(USO)
      CALL bounds_exch(VSE)

#endif  /*OBCS_TST_NEST*/
      END SUBROUTINE OBCS_TST_UVO


      SUBROUTINE OBCS_TST_UVK(UUU0,VVV0,i_step)
#ifdef OBCS_TST_NEST
! This is the TST obc for BC current.
      REAL UUU0(1:IE-I_start+1,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,1:JE-J_start+1,KE),VVV(IE,J_start:JE,KE)
      REAL M3obc_in, M3obc_out
      INTEGER inflow, i_step, obcfac, M3nudgcof

      obcfac = 4
      M3nudgcof = 1.0*86400
      M3obc_out = M3nudgcof
      M3obc_in = obcfac * M3obc_out
      
      tau = 0.0

      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)
! -----------------------  U-bcs for 4 directions ---------------------------
!--------------------------------------------------------------------------------
! Western edge
!--------------------------------------------------------------------------------
      if (WEST .and. have_g_is .and. TSTWest) then
       uko_wb0(1:3,:,:) = UUU(1:3,:,:)
       if (I_STEP .eq. 1 ) then
         ! uko_wb1 = uko_wb0
         ! obc_uNw = obcf_uw
         ! obc_uTNw = USOT(1:3,:)
       end if
       do k=1,ke  ! k
       	 gradu_wb(:,:)=0.0
         do j=je1,1,-1
            if (AMSUO(0,J,k)  .eq. 1) THEN
           ! dUdy
           ! j = je  gradu_wb = 0
           ! j = 0   gradu_wb = 0  apply !
            gradu_wb(0,j) = uko_wb1(1,j,k) - obc_uTNw(1,j) -  obc_uNw(j,k)- &
            &           uko_wb1(0,j+1,k) + obc_uTNw(0,j+1) + obc_uNw(j+1,k)
            gradu_wb(1,j) = uko_wb1(2,j,k) - obc_uTNw(2,j) -  obc_uNw(j,k)- &
            &           uko_wb1(2,j+1,k)  + obc_uTNw(2,j+1) + obc_uNw(j+1,k)
            end if
           enddo
  ! Attenetion:  the corners excluded
            do j=je1,2,-1
               if (AMSUO(0,J,k)  .eq. 1) THEN
               dUdt=uko_wb1(2,j,k) - obc_uTNw(2,j) - obc_uNw(j,k) -     &
            &           uko_wb0(2,j,k) + USOT(2,j) +  obcf_uw(j,k)
               dUdx=uko_wb0(2,j,k) - USOT(2,j) - obcf_uw(j,k) -            &
            &           uko_wb0(3,j,k) + USOT(3,j) +  obcf_uw(j,k)  ! upstream shceme

               if ((dUdt*dUdx) .lt. 0.) dUdt=0.
               if ((dUdt*(gradu_wb(1,j)+gradu_wb(1,j-1))) .gt. 0.) then  ! upstream shceme
                  dUde = gradu_wb(1,j)
               else
                  dUde = gradu_wb(1,j-1)
               endif
               cff = max(dUdx*dUdx+dUde*dUde,obcerr)
               Cx = dUdt*dUdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
                obc_uw(j,k) = (cff*(uko_wb1(1,j,k) - obc_uTNw(1,j) - obc_uNw(j,k)) + &
     &                        Cx *(uko_wb0(2,j,k) - USOT(2,j) - obcf_uw(j,k)) -        &
     &                        MAX(Ce,0.)*gradu_wb(0,j) -                &
     &                        MIN(Ce,0.)*gradu_wb(0,j-1))/               &
     &                        (cff+Cx)

               if ((dUdt*dUdx) .lt. 0.0) then
                  inflow = 1
                  tau = M3obc_in
               else
                  inflow = 0
                  tau = M3obc_out
               endif

#ifdef IMPLICIT_NUDGING
               if (tau .gt. 0.0) tau = 1.0/tau
#else
               tau = tau*dt
#endif

           if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_uw(j,k) = uko_wb1(1,j,k) +                               &
     &                   phi*(obc_uNw(j,k) - uko_wb1(1,j,k) + USOT(1,j))
#else
               obc_uw(j,k) = uko_wb1(1,j,k) +                               &
     &                   tau*(obc_uNw(j,k) - uko_wb1(1,j,k) + USOT(1,j))
#endif
           else
               obc_uw(j,k) = obc_uw(j,k)  + obcf_uw(j,k) + USOT(1,j)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
               phi = dt / (tau + dt)
               obc_uw(j,k)  = (1.0-phi)*obc_uw(j,k)  +   &
      &                               phi*obcf_uw(j,k)
#else
               obc_uw(j,k) = obc_uw(j,k) + tau*(obcf_uw(j,k) - uko_wb1(1,j,k))
#endif
#endif /*RADIATION_NUDGING*/
            endif
             obc_uNw(j,k)	 = obcf_uw(j,k)
               endif
            enddo
           obc_uNw(1,k) = obcf_uw(1,k)
           obc_uNw(je,k)= obcf_uw(je,k)
       enddo  ! k
       uko_wb1 = uko_wb0
       obc_uTNw = USOT(1:3,:)
      endif
!--------------------------------------------------------------------------------
! Eastern Edge -- U_bc
!--------------------------------------------------------------------------------
      if (EAST .and. have_g_ie .and. TSTEast) then
         uko_eb0(ie-3:ie1,:,:)=UUU(ie-3:ie1,:,:)
         if (I_STEP .eq. 1 ) then
         ! uko_eb1 = uko_eb0
         ! obc_uNe = obcf_ue
         ! obc_uTNe = USOT(ie-3:ie1,:)
         end if
         do k=1,ke,1
         	  gradu_eb(:,:)=0.0
            do j=je1,1,-1
               if (AMSUO(IE,J,K) .eq. 1) THEN
               ! dUdy
               gradu_eb(ie1,j) = uko_eb1(ie2,j,k) - obc_uTNe(ie2,j) - obc_uNe(j,k) - &
    &                           uko_eb1(ie2,j+1,k) + obc_uTNe(ie2,j+1) + obc_uNe(j+1,k)
               gradu_eb(ie, j) = uko_eb1(ie1, j,k) - obc_uTNe(ie1,j) - obc_uNe(j,k) - &
    &                           uko_eb1(ie1 ,j+1,k)  + obc_uTNe(ie1,j+1) + obc_uNe(j+1,k)
               endif
            enddo
! Attentions: the corners excluded
            do j=je2, 2,-1
               if (AMSUE(IE,J,K) .eq. 1) THEN
                dUdt= uko_eb1(ie2,j,k) - obc_uTNe(ie2,j) - obc_uNe(j,k) -        &
     &               uko_eb0(ie2,j,k) + USOT(ie2,j) +  obcf_ue(j,k)
                dUdx= uko_eb0(ie2,j,k) - USOT(ie2,j) -  obcf_ue(j,k) -       &
     &               uko_eb0(ie-3,j,k) + USOT(ie-3,j) +  obcf_ue(j,k)
                if ((dUdt*dUdx) .lt. 0.) dUdt=0.
                if ((dUdt*(gradu_eb(ie1,j)+gradu_eb(ie1,j-1))) .gt. 0.) then  ! upstream shceme
                     dUde = gradu_eb(ie1,j)
                else
                     dUde = gradu_eb(ie1,j-1)
                endif
                cff = max(dUdx*dUdx+dUde*dUde,obcerr)
                Cx = dUdt*dUdx
#ifdef RADIATION_2D
                Ce = MIN(cff,MAX(dUdt*dUde,-cff))
#else
                Ce = 0.
#endif /*RADIATION_2D*/
                obc_ue(j,k) = (cff*(uko_eb1(ie1, j,k) - obc_uTNe(ie1,j) - obc_uNe(j,k)) +   &
      &                     Cx *(uko_eb0(ie2,j,k) - USOT(ie2,j) - obc_uNe(j,k)) -          &
      &                     MAX(Ce,0.)*gradu_eb(ie,j  )-                  &
      &                     MIN(Ce,0.)*gradu_eb(ie,j-1))/                 &
      &                     (cff+Cx)

                if ((dUdt*dUdx) .lt. 0.0) then
                  inflow = 1
                  tau = M3obc_in
                else
                  tau = M3obc_out
                  inflow = 0
                endif

#ifdef IMPLICIT_NUDGING
                 IF (tau .gt. 0.0) tau=1.0/tau
#else
                 tau = tau*dt
#endif
                 if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                    phi = dt / (tau + dt)
                    obc_ue(j,k) = uko_eb1(ie,j,k) +                              &
     &                 phi*(obc_uNe(j,k) - uko_eb1(ie,j,k) + USOT(ie,j))
#else
                    obc_ue(j,k) = uko_wb1(ie,j,k) +                              &
     &                 tau*(obc_uNe(j,k) - uko_eb1(ie,j,k) + USOT(ie,j))
#endif
                 else
                    obc_ue(j,k) = obc_ue(j,k)  + obcf_ue(j,k) + USOT(ie,j)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                    phi = dt / (tau + dt)
                    obc_ue(j,k) = (1.0-phi) * obc_ue(j,k) +                      &
      &                   phi*obcf_ue(j,k)
#else
                   obc_ue(j,k) = obc_ue(j,k) + tau*(obcf_ue(j,k) - uko_eb1(ie,j,k))
#endif
#endif /*RADIATION_NUDGING*/
                  endif
                  obc_uNe(j,k) = obcf_ue(j,k)
               endif
           enddo
         obc_uNe(1,k) = obcf_ue(1,k)
         obc_uNe(je,k)= obcf_ue(je,k)
       enddo
       uko_eb1=uko_eb0
       obc_uTNe = USOT(ie-3:ie1,:)
      endif
!--------------------------------------------------------------------------------
! Southern Edge
!--------------------------------------------------------------------------------
      if (TanSouth) then
      if (SOUTH .and. have_g_je .and. TSTSouth) then
         uko_sb0=UUU(:,je2:je,:)
         if (I_STEP .eq. 1 ) then
          !  uko_sb1 = uko_sb0
          !  obc_uNs = obcf_us
          !  obc_uTNs = USOT(:,je2:je)
         endif
         do k=1,ke
         	 gradu_sb(:,:)=0.0
            do i=I_start,ie1
               if (AMSUO(I,JE,K) .eq. 1) THEN
               	! dUdx
               gradu_sb(i,je) =uko_sb1(i+1,je ,k) - obc_uTNs(i+1,je) - obc_uNs(i+1,k) -  &
     &                       uko_sb1(i,je,k) + obc_uTNs(i,je) + obc_uNs(i,k)
               gradu_sb(i,je1)=uko_sb1(i+1,je1,k) - obc_uTNs(i+1,je1) - obc_uNs(i+1,k) -  &
     &                       uko_sb1(i,je1,k) + obc_uTNs(i,je1) + obc_uNs(i,k)
              endif
            enddo
! Attentions: the corners excluded
               do i=I_start+2,ie2
                  if (AMSUO(I,JE,K) .eq. 1) THEN
                  dUdt=uko_sb1(i,je1,k) - obc_uTNs(i,je1) - obc_uNs(i,k) -        &
     &                  uko_sb0(i,je1,k) + USOT(i,je1) + obcf_us(i,k)
                  dUde=uko_sb0(i,je1,k) - USOT(i,je1) - obc_uNs(i,k) -        &
     &                  uko_sb0(i,je2,k) + USOT(i,je2) +  obcf_us(i,k)  ! upstream shceme

                  if ((dUdt*dUde) .lt. 0.) dUdt=0.
                  if ((dUdt*(gradu_sb(i-1,je1)+gradu_sb(i,je1))) .gt. 0.)then ! upstream shceme
                     dUdx = gradu_sb(i-1,je1)
                  else
                     dUdx = gradu_sb(i  ,je1)
                  end if
                  cff=MAX(dUdx*dUdx+dUde*dUde,obcerr)
                  Ce=dUdt*dUde

#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dUdt*dUdx,-cff))
#else
                  Cx = 0.
#endif
                  obc_us(i,k) = (cff*(uko_sb1(i,je, k) - obc_uTNs(i,je) - obc_uNs(i,k)) +      &
     &                           Ce *(uko_sb0(i,je1,k) - USOT(i,je1) - obcf_us(i,k)) -         &
     &                           MAX(Cx,0.)*gradu_sb(i-1,je)-            &
     &                           MIN(Cx,0.)*gradu_sb(i  ,je))/           &
     &                           (cff+Ce)

                 if ((dUdt*dUde) .lt. 0.0) then
                     inflow = 1
                     tau = M3obc_in
                 else
                     tau = M3obc_out
                     inflow = 0
                 end if

#ifdef IMPLICIT_NUDGING
                 if (tau.gt.0.0) tau = 1.0/tau
#else
                 tau = tau*dt
#endif

                  if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                   phi = dt / (tau + dt)
                   obc_us(i,k) = uko_sb1(i,je,k) +                              &
     &                 phi*(obc_uNs(i,k) - uko_sb1(i,je,k) + USOT(i,je))
#else
                  obc_us(i,k) = uko_sb1(i,je,k) +                              &
     &                 tau*(obc_uNs(i,k) - uko_sb1(i,je,k) + USOT(i,je))
#endif
                 else
                  obc_us(i,k) = obc_us(i,k)  + obcf_us(i,k) + USOT(i,je)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                  phi = dt / (tau + dt)
                  obc_us(i,k) = (1.0-phi) * obc_us(i,k) +                      &
     &                   phi*obcf_us(i,k)
#else
                  obc_us(i,k) = obc_us(i,k) + tau*(obcf_us(i,k) - uko_sb1(i,je,k))
#endif
#endif /*RADIATION_NUDGING*/
                  endif
                 obc_uNs(i,k) = obcf_us(i,k)
                 endif
               enddo
               obc_uNs(I_start:I_start+1,k) = obcf_us(I_start:I_start+1,k)
               obc_uNs(ie1:ie,k)= obcf_us(ie1:ie,k)
         enddo
         uko_sb1=uko_sb0
         obc_uTNs = USOT(:,je2:je)
      endif
      endif
!--------------------------------------------------------------------------------
! Northern Edge
!--------------------------------------------------------------------------------
      if (TanNorth) then
      if (NORTH .and. have_g_js .and. TSTNorth) then
         uko_nb0=UUU(:,1:3,:)
         if (I_STEP .eq. 1 ) then
           ! uko_nb1 = uko_nb0
           ! obc_uNn = obcf_un
           ! obc_uTNn = USOT(:,1:3)
         endif
         do k=1,ke
         gradu_nb(:,:)=0.0
            do i=I_start,ie1
               if (AMSUO(I,1,K) .eq. 1) THEN
               	! dUdx
               gradu_nb(i,2) = uko_nb1(i+1,2,k)  - obc_uTNn(i+1,2) - obc_uNn(i+1,k) -  &
     &                      uko_nb1(i,2,k) + obc_uTNn(i,2) + obc_uNn(i,k)
               gradu_nb(i,1) = uko_nb1(i+1,1,k) - obc_uTNn(i+1,1) -  obc_uNn(i+1,k) -  &
     &                      uko_nb1(i,1,k) + obc_uTNn(i,1) + obc_uNn(i,k)
               ENDIF
             ENDDO

! Attentions: the corners excluded
              do i=I_start+2,ie2
                  if (AMSUO(I,1,K) .eq. 1) THEN
                  dUdt=uko_nb1(i,2,k) - obc_uTNn(i,2) -  obc_uNn(i,k) -           &
     &                      uko_nb0(i,2,k) + USOT(i,2) +  obcf_un(i,k)
                  dUde=uko_nb0(i,2,k) - USOT(i,2) -  obcf_un(i,k) -           &
     &                      uko_nb0(i,3,k) + USOT(i,3) +  obcf_un(i,k)  ! upstream shceme

                  if ((dUdt*dUde) .lt. 0.) dUdt=0.
                  if ((dUdt*(gradu_nb(i-1,2)+gradu_nb(i,2))) .gt. 0.) then
                     dUdx = gradu_nb(i-1,2)
                  else
                     dUdx = gradu_nb(i  ,2)
                  end if
                  cff=MAX(dUdx*dUdx+dUde*dUde,obcerr)
                  Ce=dUdt*dUde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dUdt*dUdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                 obc_un(i,k) = (cff*(uko_nb1(i,1,k) - obc_uTNn(i,1) -obc_uNn(i,k)) +     &
     &                          Ce *(uko_nb0(i,2,k) - USOT(i,2) - obcf_un(i,k))  -          &
     &                          MAX(Cx,0.)*gradu_nb(i-1,1)-              &
     &                          MIN(Cx,0.)*gradu_nb(i  ,1))/             &
     &                          (cff+Ce)

                 if ((dUdt*dUde) .lt. 0.0) then
                   inflow = 1
                   tau = M3obc_in
                else
                   tau = M3obc_out
                   inflow = 0
                end if

#ifdef IMPLICIT_NUDGING
                IF (tau .gt. 0.0) tau=1.0/tau
#else
                tau = tau*dt
#endif
                if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                	phi = dt / (tau + dt)
                  obc_un(i,k) = uko_nb1(i,1,k) +                               &
     &                 phi*(obc_uNn(i,k) - uko_nb1(i,1,k) + USOT(i,1))
#else
                   obc_un(i,k) = uko_nb1(i,1,k) +                               &
     &                 tau*(obc_uNn(i,k) - uko_nb1(i,1,k) + USOT(i,1))
#endif
                else
                    obc_un(i,k) = obc_un(i,k)  + obcf_un(i,k) + USOT(i,1)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                     phi = dt / (tau + dt)
                    obc_un(i,k) = (1.0-phi) * obc_un(i,k) +                      &
      &                   phi*obcf_un(i,k)
#else
                    obc_un(i,k) = obc_un(i,k) + tau*(obcf_un(i,k) - uko_nb1(i,1,k))
#endif
#endif /*RADIATION_NUDGING*/
                endif
                obc_uNn(i,k) = obcf_un(i,k)
              endif
            end do ! i
         obc_uNn(I_start:I_start+1,k) = obcf_un(I_start:I_start+1,k)
         obc_uNn(ie1:ie,k)  = obcf_un(ie1:ie,k)
       end do ! k
       uko_nb1=uko_nb0
       obc_uTNn = USOT(:,1:3)
      endif
      endif

! -----------------------  V-bcs for 4 directions ---------------------------
!-----------------------------------------------------------------------------
! Southern Edge -- V_bc
!-----------------------------------------------------------------------------

      if (SOUTH .and. have_g_je .and. TSTSouth) then
         vke_sb0(:,je-3:je-1,:)=VVV(:,je-3:je-1,:)
         if (I_STEP .eq. 1 ) then
         ! vke_sb1 = vke_sb0
         ! obc_vNs = obcf_vs
         ! obc_vTNs = VSET(:,je-3:je-1)
         end if
         do k=1,ke
            gradv_sb(:,:)=0.0
            do i=2,ie
               if (AMSUE(I,JE,K) .eq. 1) THEN
               	! dVdx
               gradv_sb(i,je)   = vke_sb1(i,je1, k) - obc_vTNs(i,je1) - obc_vNs(i,k) -    &
     &                     vke_sb1(i-1,je1, k) + obc_vTNs(i-1,je1)+obc_vNs(i-1,k)
               gradv_sb(i,je1)  = vke_sb1(i,je2,k) - obc_vTNs(i,je2) - obc_vNs(i,k) -    &
     &                     vke_sb1(i-1,je2,k) + obc_vTNs(i-1,je2)+obc_vNs(i-1,k)
               endif
            enddo
! Attentions: the corners excluded
            do i=2,ie1
                if (AMSUE(I,JE,K) .eq. 1) THEN
                  dVdt=vke_sb1(i,je2,k) - obc_vTNs(i,je2) - obc_vNs(i,k) -    &
     &                     vke_sb0(i,je2,k) + VSET(i,je2) + obcf_vs(i,k)
                  dVde=vke_sb0(i,je2,k) - VSET(i,je2) - obc_vNs(i,k) -    &
     &                     vke_sb0(i,je-3,k) + VSET(i,je-3) + obcf_vs(i,k)   ! upstream shceme
               if ((dVdt*dVde).lt.0.) dVdt=0.
               if ((dVdt*(gradv_sb(i,je1)+gradv_sb(i+1,je1))).gt.0.) then ! upstream shceme
                  dVdx=gradv_sb(i,je1)
               else
                  dVdx=gradv_sb(i+1,je1)
               end if
               cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
               Ce=dVdt*dVde
#ifdef RADIATION_2D
               Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#else
               Cx=0.
#endif /*RADIATION_2D*/
               obc_vs(i,k) = (cff*(vke_sb1(i,je1,k) - obc_vTNs(i,je1) - obc_vNs(i,k) )+                      &
     &                        Ce *(vke_sb0(i,je2,k)  - VSET(i,je2) - obcf_vs(i,k))-                     &
     &                        MAX(Cx,0.)*gradv_sb(i,je)-                 &
     &                        MIN(Cx,0.)*gradv_sb(i+1,je))/              &
     &                        (cff+Ce)
               if ((dVdt*dVde) .lt. 0.0) then
                  inflow = 1
                  tau = M3obc_in
               else
                  tau = M3obc_out
                  inflow = 0
               end if

#ifdef IMPLICIT_NUDGING
               IF (tau .gt. 0.0) tau=1.0/tau
#else
                tau = tau*dt
#endif
                if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                    phi = dt / (tau + dt)
                    obc_vs(i,k) = vke_sb1(i,je1,k) +                              &
       &                 phi*(obc_vNs(i,k) - vke_sb1(i,je1,k) + VSET(i,je1))
#else
                    obc_vs(i,k) = vke_sb1(i,je1,k) +                              &
       &                 tau*(obc_vNs(i,k) - vke_sb1(i,je1,k) + VSET(i,je1))
#endif
                else
                    obc_vs(i,k) = obc_vs(i,k)  + obcf_vs(i,k) + VSET(i,je1)

#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                    phi = dt / (tau + dt)
                   obc_vs(i,k) = (1.0-phi) * obc_vs(i,k) +                      &
      &                   phi*obcf_vs(i,k)
#else
                 obc_vs(i,k) = obc_vs(i,k) + tau*(obcf_vs(i,k) - vke_sb1(i,je1,k))
#endif
#endif /*RADIATION_NUDGING*/
                endif
                obc_vNs(i,k) = obcf_vs(i,k)
              endif
            enddo
            obc_vNs(1,k) = obcf_vs(1,k)
           obc_vNs(ie,k)= obcf_vs(ie,k)
          enddo
          vke_sb1=vke_sb0
          obc_vTNs = VSET(:,je-3:je1)
      endif
!---------------------------------------------------------------------------------
! Northern edge.
!---------------------------------------------------------------------------------
      if (NORTH .and. have_g_js .and. TSTNorth) then
         ! vke_nb0=VVV(:,0:2,:)
         vke_nb0(:,1:3,:)=VVV(:,1:3,:)
         if (I_STEP .eq. 1 ) then
         ! vke_nb1 = vke_nb0
         ! obc_vNn = obcf_vn
         ! obc_vTNn = VSET(:,1:3)
         end if
         do k=1,ke
            gradv_nb(:,:)=0.0
            do i=2,ie
               if (AMSUE(I,0,K) .eq. 1) THEN
               gradv_nb(i,1) = vke_nb1(i,2,k)- obc_vTNn(i,2)-obc_vNn(i,k) - &
     &                      vke_nb1(i-1,2,k) + obc_vTNn(i-1,2) + obc_vNn(i,k)
               gradv_nb(i,0) = vke_nb1(i,1,k) -  obc_vTNn(i,1) - obc_vNn(i,k) -  &
     &                      vke_nb1(i-1,1,k) + obc_vTNn(i,1) + obc_vNn(i,k)
               end if
            end do
! Attentions: the corners excluded
            do i=2,ie1
                  if (AMSUE(I,0,K) .eq. 1) THEN
                    dVdt=vke_nb1(i,2,k)  - obc_vTNn(i,2) -  obc_vNn(i,k) -    &
     &                      vke_nb0(i,2,k) + VSET(i,2) +  obcf_vn(i,k)
                    dVde=vke_nb0(i,2,k)  - VSET(i,2) -  obcf_vn(i,k) -           &
     &                      vke_nb0(i,3,k) + VSET(i,3) +  obcf_vn(i,k)  ! upstream shceme
                  if ((dVdt*dVde) .lt. 0.) dVdt=0.
                  if ((dVdt*(grad_nb(i,1)+grad_nb(i+1,1))) .gt. 0.) then
                       dVdx = grad_nb(i  ,1)
                 else
                       dVdx = grad_nb(i+1,1)
                 endif
                 cff = max(dVdx*dVdx+dVde*dVde,obcerr)
                 Ce = dVdt*dVde
#ifdef RADIATION_2D
                 Cx = MIN(cff,MAX(dVdt*dVdx,-cff))
#else
                 Cx = 0.
#endif /*RADIATION_2D*/
                  obc_vn(i,k) = (cff*(vke_nb1(i,1,k)  -obc_vTNn(i,1) - obc_vNn(i,k)) + &
      &                           Ce *(vke_nb0(i,2,k)  - VSET(i,2) - obcf_vn(i,k)) -             &
      &                           MAX(Cx,0.)*gradv_nb(i,  0)-             &
      &                           MIN(Cx,0.)*gradv_nb(i+1,0))/            &
      &                           (cff+Ce)
                  if ((dVdt*dVde) .lt. 0.0) then
                      inflow = 1
                       tau = M3obc_in
                  else
                       tau = M3obc_out
                       inflow = 0
                  endif
#ifdef IMPLICIT_NUDGING
                  IF (tau .gt. 0.0) tau=1.0/tau
#else
                  tau = tau*dt
#endif
                   if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                     phi = dt / (tau + dt)
                     obc_vn(i,k) = vke_nb1(i,1,k) +                               &
         &                 phi*(obc_vNn(i,k) - vke_nb1(i,1,k) + VSET(i,1))
#else
                     obc_vn(i,k) = vke_nb1(i,1,k) +                               &
        &                 tau*(obc_vNn(i,k) - vke_nb1(i,1,k) + VSET(i,1))
#endif
                   else
                     obc_vn(i,k) = obc_vn(i,k)  + obcf_vn(i,k) + VSET(i,1)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                     phi = dt / (tau + dt)
                     obc_vn(i,k) = (1.0-phi) * obc_vn(i,k) +                    &
       &                   phi*obcf_vn(i,k)
#else
                     obc_vn(i,k) = obc_vn(i,k) + tau*(obcf_vn(i,k) - vke_nb1(i,1,k))
#endif
#endif /*RADIATION_NUDGING*/
                   endif
                   obc_vNn(i,k) = obcf_vn(i,k)
                  endif
           enddo ! i
         obc_vNn(1,k) = obcf_vn(1,k)
         obc_vNn(ie,k) = obcf_vn(ie,k)
        enddo ! k
        vke_nb1=vke_nb0
        obc_vTNn = VSET(:,1:3)
      endif

!--------------------------------------------------------------------------------
!  Western edge
!--------------------------------------------------------------------------------
      if (TanWest) then
      if (WEST .and. have_g_is .and. TSTWest) then
         vke_wb0=VVV(1:3,:,:)
         if (I_STEP .eq. 1 ) then
          ! vke_wb1 = vke_wb0
          ! obc_vNw = obcf_vw
          ! obc_vTNw = VSET(1:3,:)
         endif
         do k=1,ke
            gradv_wb(:,:)=0.0
            do j=je,J_start+1,-1 ! ATTEN: The direction is from S -> N.
               if (AMSUE(1,J,K) .eq. 1) THEN
               	! dVdy
               gradv_wb(1,j) =  vke_wb1(1,j-1,k) - obc_vTNw(1,j-1) - obc_vNw(j-1,k) - &
     &                     - vke_wb1(1,j,k) +  obc_vTNw(1,j) + obc_vNw(j,k)
               gradv_wb(2,j) =  vke_wb1(2,j-1,k) - obc_vTNw(2,j-1) - obc_vNw(j-1,k) - &
     &                     - vke_wb1(2,j,k) +  obc_vTNw(2,j) + obc_vNw(j,k)
               end if
             enddo
! Attentions: the corners excluded
             do j=je2,J_start+2,-1
                  if (AMSUE(1,J,K) .eq. 1) THEN

                  dVdt=vke_wb1(2,j,k) - obc_vTNw(2,j) - obc_vNw(j,k) -   &
     &               vke_wb0(2,j,k)  + VSET(2,j) +  obcf_vw(j,k)
                  dVdx=vke_wb0(2,j,k) - VSET(2,j) -  obcf_vw(j,k) -           &
     &                vke_wb0(3,j,k)  + VSET(3,j) +  obcf_vw(j,k)  ! upstream shceme

                  if ((dVdt*dVdx).lt.0.) dVdt=0.
                  if ((dVdt*(gradv_wb(2,j+1)+gradv_wb(2,j))).gt.0.) then ! upstream shceme
                     dVde=gradv_wb(2,j+1)
                  else
                     dVde=gradv_wb(2,j  )
                  end if
                     cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
                     Cx=dVdt*dVdx
#ifdef RADIATION_2D
                     Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
                     Ce=0.
#endif /*RADIATION_2D*/
                     obc_vw(j,k) = (cff*(vke_wb1(1,j,k) - obc_vTNw(1,j) - obc_vNw(j,k)  ) +     &
     &                             Cx *(vke_wb0(2,j,k) -  VSET(2,j) - obcf_vw(j,k) )-                  &
     &                             MAX(Ce,0.)*gradv_wb(1,j+1)-           &
     &                             MIN(Ce,0.)*gradv_wb(1,j  ))/          &
     &                             (cff+Cx)
                      if ((dVdt*dVdx) .lt. 0.0) then
                        inflow = 1
                        tau = M3obc_in
                      else
                        tau = M3obc_out
                        inflow = 0
                      endif
#ifdef IMPLICIT_NUDGING
                      IF (tau .gt. 0.0) tau = 1.0/tau
#else
                      tau = tau*dt
#endif
                      if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                      	phi = dt / (tau + dt)
                      	obc_vw(j,k) = vke_wb1(1,j,k) +                               &
     &                   phi*(obc_vNw(j,k) - vke_wb1(1,j,k) + VSET(1,j))
#else
                        obc_vw(j,k) = vke_wb1(1,j,k) +                               &
     &                   tau*(obc_vNw(j,k) - vke_wb1(1,j,k) + VSET(1,j))
#endif
                      else
                      	 obc_vw(j,k) = obc_vw(j,k)  + obcf_vw(j,k) + VSET(1,j)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                       	 phi = dt / (tau + dt)
                         obc_vw(j,k) = (1.0-phi) * obc_vw(j,k) +   &
      &                   phi*obcf_vw(j,k)
#else
                         obc_vw(j,k) = obc_vw(j,k) + tau*(obcf_vw(j,k) - vke_wb1(1,j,k))
#endif
#endif /*RADIATION_NUDGING*/
                      endif
                      obc_vNw(j,k) = obcf_vw(j,k)
                  endif
             enddo
             obc_vNw(J_start:J_start+1,k) = obcf_vw(J_start:J_start+1,k)
             obc_vNw(je1:je,k)= obcf_vw(je1:je,k)
        enddo
        vke_wb1 = vke_wb0
        obc_vTNw = VSET(1:3,:)
      endif
      endif

!-------------------------------------------------------------------------
!  Eastern edge
!-------------------------------------------------------------------------
      if (TanEast) then
      if (EAST .and. have_g_ie .and. TSTEast) then
         vke_eb0=VVV(ie2:ie,:,:)
         if (I_STEP .eq. 1 ) then
          ! vke_eb1 = vke_eb0
          ! obc_vNe = obcf_ve
          ! obc_vTNe = VSET(ie2:ie,:)
         end if
         do k=1,ke
            gradv_eb(:,:)=0.0
            do j=je,J_start+1,-1
               if (AMSUE(IE,J,K) .eq. 1) THEN
               	! dVdy
                  gradv_eb(ie1,j) = vke_eb1(ie1,j-1,k) - obc_vTNe(ie1,j-1) - obc_vNe(j-1,k) -  &
     &                       vke_eb1(ie1,j,k) + obc_vTNe(ie1,j)  + obc_vNe(j,k)
                  gradv_eb(ie, j) = vke_eb1(ie, j-1,k)  - obc_vTNe(ie,j-1) - obc_vNe(j-1,k) -  &
     &                       vke_eb1(ie, j,k) + obc_vTNe(ie,j)+obc_vNe(j,k)
               endif
             enddo
! Attentions: the corners excluded
             do j=je2,J_start+2,-1
                  if (AMSUE(IE,J,K) .eq. 1) THEN
                  dVdt=vke_eb1(ie1,j,k) - obc_vTNe(ie1,j-1) -  obc_vNe(j-1,k) -          &
     &                       vke_eb0(ie1,j,k) + VSET(ie1,j) + obc_vNe(j,k)
                  dVdx=vke_eb0(ie1,j,k) - VSET(ie,j-1) -  obc_vNe(j-1,k) -          &
     &                       vke_eb0(ie2,j,k) + VSET(ie,j) + obc_vNe(j,k)

                  if ((dVdt*dVdx).lt.0.) dVdt=0.
                  if ((dVdt*(gradv_eb(ie1,j+1)+gradv_eb(ie1,j))).gt.0.)then
                     dVde=gradv_eb(ie1,j+1)
                  else
                     dVde=gradv_eb(ie1,j  )
                  end if
                     cff=MAX(dVdx*dVdx+dVde*dVde,obcerr)
                     Cx=dVdt*dVdx
#ifdef RADIATION_2D
                     Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
                     Ce=0.
#endif /*RADIATION_2D*/
                     obc_ve(j,k)  = (cff*(vke_eb1(ie,j,k) - obc_vTNe(ie,j) - obc_vNe(j,k)  )+   &
     &                               Cx *(vke_eb0(ie1,j,k) - VSET(ie1,j) - obcf_ve(j,k)) -       &
     &                               MAX(Ce,0.)*gradv_eb(ie,j+1)-        &
     &                               MIN(Ce,0.)*gradv_eb(ie,j  ))/       &
     &                               (cff+Cx)

                     if ((dVdt*dVdx) .lt. 0.0) then
                      inflow = 1
                      tau = M3obc_in
                     else
                      tau = M3obc_out
                      inflow = 0
                     endif
#ifdef IMPLICIT_NUDGING
                      IF (tau .gt. 0.0) tau = 1.0/tau
#else
                      tau = tau*dt
#endif
                      if (inflow .eq. 1) then
#ifdef IMPLICIT_NUDGING
                        phi = dt / (tau + dt)
                        obc_ve(j,k) = vke_eb1(ie,j,k) +     &
     &                     phi*(obc_vNe(j,k) - vke_eb1(ie,j,k) + VSET(ie,j))
#else
                        obc_ve(j,k) = vke_wb1(ie,j,k) +    &
     &                      tau*(obc_vNe(j,k) - vke_eb1(ie,j,k) + VSET(ie,j))
#endif
                      else
                         obc_ve(j,k) = obc_ve(j,k)  + obcf_ve(j,k) + VSET(ie,j)
#ifdef RADIATION_NUDGING
#ifdef IMPLICIT_NUDGING
                          phi = dt / (tau + dt)
                          obc_ve(j,k) = (1.0-phi) * obc_ve(j,k) +    &
      &                       phi*obcf_ve(j,k)
#else
                          obc_ve(j,k) = obc_ve(j,k) + tau*(obcf_ve(j,k) - vke_eb1(ie,j,k))
#endif
#endif /*RADIATION_NUDGING*/
                       endif
                  obc_vNe(j,k) = obcf_ve(j,k)
                  endif
             enddo
             obc_vNe(J_start:J_start+1,k) = obcf_ve(J_start:J_start+1,k)
             obc_vNe(je1:je,k)= obcf_ve(je1:je,k)
         enddo
         vke_eb1 = vke_eb0
         obc_vTNe = VSET(ie2:ie,:)
      endif
      endif

! Update OBC predict uv to model fields

      CALL UPDATE_ORLANSKI_UV(UUU,VVV) ! in mo_obcs.f90

      UUU0(:,:,:)=UUU(:,:,:)
      VVV0(:,:,:)=VVV(:,:,:)
#endif /*OBCS_TST_NEST*/
      END SUBROUTINE OBCS_TST_UVK

      SUBROUTINE OBCS_TST_TS(TTT0,SSS0,i_step)
#ifdef OBCS_TST_NEST
      REAL TTT(IE,JE,KE),TTT0(IE,JE,KE),SSS(IE,JE,KE),SSS0(IE,JE,KE)
      REAL tobc_in, tobc_out, tnudgcof
      INTEGER obcfac
      phi = 0.0
      tau = 0.0
      obcfac = 4.0
      tnudgcof = 1.0*86400
      
      tobc_out = tnudgcof
      tobc_in   = obcfac * tobc_out

      TTT(:,:,:)=TTT0(:,:,:)
      SSS(:,:,:)=SSS0(:,:,:)
! Western edge
      if (WEST .and. have_g_is .and. TSTWest ) then
       tho_wb0 = TTT(1:3,:,:)
       do k=1,ke
         grad_wb=0.0
         do j=je1,1,-1 ! ATTEN: The mpiom's je is from north to south,oppsite to Roms.
         	 ! dTdy
            grad_wb(1,j) = tho_wb1(1,j,k) - cli_tNw(j,k) -   &
     &                     tho_wb1(1,j+1,k)+ cli_tNw(j+1,k)
            grad_wb(2,j) = tho_wb1(2,j,k)  - cli_tNw(j,k) -  &
     &                     tho_wb1(2,j+1,k) + cli_tNw(j+1,k)
         enddo
! Attentions: the corners excluded
            do j=je1,2,-1
            	if  (weto(1,j,k) .eq. 1) then
               dTdt=tho_wb1(2,j,k) - cli_tNw(j,k) -                     &
     &              tho_wb0(2,j,k) + cli_tw(j,k)
               dTdx=tho_wb0(2,j,k) - cli_tw(j,k)  -                     &
     &              tho_wb0(3,j,k)  + cli_tw(j,k)  ! upstream scheme

               if ((dTdt*dTdx) .lt. 0.) dTdt=0.
               if ((dTdt*(grad_wb(2,j)+grad_wb(2,j-1))) .gt. 0.) then
                  dTde = grad_wb(2,j)
               else
                  dTde = grad_wb(2,j-1)
               end if
               cff = max(dTdx*dTdx+dTde*dTde,obcerr)
               Cx = dTdt*dTdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dTdt*dTde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
               obc_tw(j,k) = (cff*(tho_wb1(1,j,k)- cli_tNw(j,k) ) +  &
     &                        Cx *(tho_wb0(2,j,k)  - cli_tw(j,k)) -      &
     &                        MAX(Ce,0.)*grad_wb(1,j  )-                &
     &                        MIN(Ce,0.)*grad_wb(1,j-1))/               &
     &                        (cff+Cx)

               if ((dTdt*dTdx) .lt. 0.0) then
                  inflow = 1
                  tau = tobc_in
               else
                  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

               if (inflow .eq. 1) then
                   obc_tw(j,k) = tho_wb1(1,j,k) + tau *                 &
     &                           (cli_tNw(j,k) - tho_wb1(1,j,k))
               else
                   obc_tw(j,k) = obc_tw(j,k) + cli_tw(j,k)
#ifdef RADIATION_NUDGING
                   obc_tw(j,k) = obc_tw(j,k) +      &
     &                       tau*(cli_tw(j,k)-             &
     &                            tho_wb1(1,j,k) )
#endif /*RADIATION_NUDGING*/
               endif
               cli_tNw(j,k) = cli_tw(j,k)
              endif
            enddo
            cli_tNw(je,k) = cli_tw(je,k)
            cli_tNw(1, k) = cli_tw(1, k)
       enddo
       tho_wb1=tho_wb0
      end if

!---------------------------------------------------------------------------
!  Eastern edge
!---------------------------------------------------------------------------
!
      if (EAST .and. have_g_ie .and. TSTEast ) then
         tho_eb0=TTT(ie2:ie,:,:)
         do k=1,ke
           grad_eb=0.0
            do j=je1,1,-1
               grad_eb(ie1,j) = tho_eb1(ie1,j,k) - cli_tNe(j,k) -       &
     &                          tho_eb1(ie1,j+1,k) + cli_tNe(j+1,k)
               grad_eb(ie, j) = tho_eb1(ie, j,k) - cli_tNe(j,k) -       &
     &                          tho_eb1(ie ,j+1,k) + cli_tNe(j+1,k)
            end do
! Attentions: the corners excluded
            do j=je1,2,-1
            	if (weto(ie,j,k) .eq. 1) then
               dTdt=tho_eb1(ie1,j,k) - cli_tNe(j,k) -       &
     &                          tho_eb0(ie1,j,k)  + cli_te(j,k)
               dTdx=tho_eb0(ie1,j,k) - cli_te(j,k) -       &
     &                          tho_eb0(ie2,j,k)  + cli_te(j,k)  ! upstream scheme
            if ((dTdt*dTdx) .lt. 0.) dTdt=0.
            if ((dTdt*(grad_eb(ie1,j)+grad_eb(ie1,j-1))) .gt. 0.) then
               dTde = grad_eb(ie1,j)
            else
               dTde = grad_eb(ie1,j-1)
            end if
            cff = max(dTdx*dTdx+dTde*dTde,obcerr)
            Cx = dTdt*dTdx
#ifdef RADIATION_2D
            Ce = MIN(cff,MAX(dTdt*dTde,-cff))
#else
            Ce = 0.
#endif /*RADIATION_2D*/
            obc_te(j,k) = (cff*(tho_eb1(ie, j,k) - cli_tNe(j,k)) +      &
     &                     Cx *(tho_eb0(ie1,j,k) - cli_te(j,k )) -      &
     &                     MAX(Ce,0.)*grad_eb(ie,j  ) -                 &
     &                     MIN(Ce,0.)*grad_eb(ie,j-1))/                 &
     &                     (cff+Cx)

               if ((dTdt*dTdx) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
               	  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

               if (inflow .eq. 1) then
                 obc_te(j,k) = tho_eb1(ie,j,k) + tau *                   &
     &                           (cli_tNe(j,k) - tho_eb1(ie,j,k))
               else
                 obc_te(j,k) = obc_te(j,k) + cli_te(j,k)
#ifdef RADIATION_NUDGING
               obc_te(j,k) = obc_te(j,k) +                              &
     &                       tau*(cli_te(j,k)-                          &
     &                            tho_eb1(ie,j,k) )

#endif /*RADIATION_NUDGING*/
                endif
                cli_tNe(j,k) = cli_te(j,k)
               endif
            end do
            cli_tNe(je,k) = cli_te(je,k)
            cli_tNe(1, k) = cli_te(1, k)
         end do
         tho_eb1=tho_eb0
      endif
!----------------------------------------------------------------------------
!  Southern edge
!----------------------------------------------------------------------------
      if (SOUTH .and. have_g_je .and. TSTSouth ) then
         tho_sb0=TTT(:,je2:je,:)
         do k=1,ke
            grad_sb=0.0
            do i=2,ie
               grad_sb(i,je1) = tho_sb1(i,je1 ,k) - cli_tNs(i,k) -      &
     &                          tho_sb1(i-1,je1 ,k) + cli_tNs(i-1,k)
               grad_sb(i,je)  = tho_sb1(i,je,k)  - cli_tNs(i,k) -        &
     &                          tho_sb1(i-1,je,k) + cli_tNs(i-1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
               	 if (weto(i,je,k) .eq. 1) then
                  dTdt=tho_sb1(i,je1,k) -  cli_tNs(i,k) -                &
     &                          tho_sb0(i,je1,k) + cli_ts(i,k)
                  dTde=tho_sb0(i,je1,k) - cli_ts(i,k) -                 &
     &                          tho_sb0(i,je2,k)  + cli_ts(i,k)  ! upstream scheme

                  if ((dTdt*dTde) .lt. 0.) dTdt=0.
                  if ((dTdt*(grad_sb(i,je1)+grad_sb(i+1,je1))) .gt. 0.)then
                     dTdx = grad_sb(i,je1)
                  else
                     dTdx = grad_sb(i+1,je1)
                  end if
                  cff=MAX(dTdx*dTdx+dTde*dTde,obcerr)
                  Ce=dTdt*dTde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dTdt*dTdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_ts(i,k) = (cff*(tho_sb1(i,je, k)-cli_tNs(i,k)) +  &
     &                           Ce *(tho_sb0(i,je1,k)-cli_ts(i,k))  -  &
     &                           MAX(Cx,0.)*grad_sb(i,je)   -           &
     &                           MIN(Cx,0.)*grad_sb(i+1,je))/           &
     &                           (cff+Ce)

               if ((dTdt*dTde) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
               	  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt
               if (inflow .eq. 1) then
                     obc_ts(i,k) = tho_sb1(i,je,k) + tau *              &
     &                             (cli_tNs(i,k) - tho_sb1(i,je,k))
               else
                     obc_ts(i,k) = obc_ts(i,k) + cli_ts(i,k)
#ifdef RADIATION_NUDGING
               obc_ts(i,k) = obc_ts(i,k) +                              &
     &                       tau*(cli_ts(i,k)-                          &
     &                            tho_sb1(i,je,k) )
#endif /*RADIATION_NUDGING*/
                endif
                cli_tNs(i,k) = cli_ts(i,k)
               endif
               end do
               cli_tNs(ie,k) = cli_ts(ie,k)
               cli_tNs(1, k) = cli_ts(1, k)
         end do
         tho_sb1=tho_sb0
      endif

!
!-----------------------------------------------------------------------------
!  Northern edge
!-----------------------------------------------------------------------------

      if (NORTH .and. have_g_js .and. TSTNorth ) then
         tho_nb0=TTT(:,1:3,:)
         do k=1,ke
            grad_nb=0.0
            do i=2,ie
               grad_nb(i,2) = tho_nb1(i,2,k) - cli_tNn(i,k) -           &
     &                        tho_nb1(i-1,2,k) + cli_tNn(i-1,k)
               grad_nb(i,1) = tho_nb1(i,1,k) - cli_tNn(i,k) -           &
     &                        tho_nb1(i-1,1,k) + cli_tNn(i-1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
               	 if( weto(i,1,k) .eq. 1) then
                  dTdt=tho_nb1(i,2,k) - cli_tNn(i,k) -                  &
     &                 tho_nb0(i,2,k)  + cli_tn(i,k)
                  dTde=tho_nb0(i,2,k) - cli_tn(i,k)  -                  &
     &                 tho_nb0(i,3,k) + cli_tn(i,k) ! upstream scheme

                  if ((dTdt*dTde) .lt. 0.) dTdt=0.
                  if ((dTdt*(grad_nb(i,2)+grad_nb(i+1,2))) .gt. 0.) then
                     dTdx = grad_nb(i,2)
                  else
                     dTdx = grad_nb(i+1,2)
                  end if
                  cff=MAX(dTdx*dTdx+dTde*dTde,obcerr)
                  Ce=dTdt*dTde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dTdt*dTdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_tn(i,k) = (cff*(tho_nb1(i,1,k)-cli_tNn(i,k)) +    &
     &                          Ce *(tho_nb0(i,2,k) -cli_tn(i,k) ) -    &
     &                          MAX(Cx,0.)*grad_nb(i  ,1)-              &
     &                          MIN(Cx,0.)*grad_nb(i+1,1))/             &
     &                          (cff+Ce)
               if ((dTdt*dTde) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
               	  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

              if (inflow .eq. 1) then
                     obc_tn(i,k) = tho_nb1(i,1,k) + tau *               &
     &                             (cli_tNn(i,k) - tho_nb1(i,1,k))
              else
                     obc_tn(i,k) = obc_tn(i,k) + cli_tn(i,k)
#ifdef RADIATION_NUDGING
                     obc_tn(i,k) = obc_tn(i,k) +                  &
     &                       tau*(cli_tn(i,k)-                          &
     &                            tho_nb1(i,1,k) )

#endif /*RADIATION_NUDGING*/
               endif
               cli_tNn(i,k) = cli_tn(i,k)
               endif
               end do
               cli_tNn(ie,k) = cli_tn(ie,k)
               cli_tNn(1, k) = cli_tn(1, k)
         end do
         tho_nb1=tho_nb0
      endif

!---------------------------------------------------------------------------------------------
!========================== Sanlity ============================!
!---------------------------------------------------------------------------------------------

! Western edge
      if (WEST .and. have_g_is .and. TSTWest) then
       sao_wb0 = SSS(1:3,:,:)
       do k=1,ke
         grad_wb=0.0
         do j=je1,1,-1 ! ATTEN: The mpiom's je is from north to south,oppsite to Roms.
            grad_wb(1,j) = sao_wb1(1,j,k) - cli_sNw(j,k) -              &
     &             sao_wb1(1,j+1,k)  + cli_sNw(j+1,k)
            grad_wb(2,j) = sao_wb1(2,j,k) -  cli_sNw(j,k) -              &
     &             sao_wb1(2,j+1,k)  + cli_sNw(j+1,k)
         enddo
! Attentions: the corners excluded
            do j=je1,2,-1
            	if( weto(1,j,k) .eq. 1) then
               dSdt=sao_wb1(2,j,k) - cli_sNw(j,k) -                   &
     &             sao_wb0(2,j,k) + cli_sw(j,k)
               dSdx=sao_wb0(2,j,k) - cli_sw(j,k) -                     &
     &             sao_wb0(3,j,k) + cli_sw(j,k)

               if ((dSdt*dSdx) .lt. 0.) dSdt=0.
               if ((dSdt*(grad_wb(2,j)+grad_wb(2,j-1))) .gt. 0.) then
                  dSde = grad_wb(2,j)
               else
                  dSde = grad_wb(2,j-1)
               end if
               cff = max(dSdx*dSdx+dSde*dSde,obcerr)
               Cx = dSdt*dSdx
#ifdef RADIATION_2D
               Ce = MIN(cff,MAX(dSdt*dSde,-cff))
#else
               Ce = 0.
#endif /*RADIATION_2D*/
               obc_sw(j,k) = (cff*(sao_wb1(1,j,k) - cli_sNw(j,k) ) +    &
     &                        Cx *(sao_wb0(2,j,k) - cli_sw(j,k)  ) -    &
     &                        MAX(Ce,0.)*grad_wb(1,j  )-                &
     &                        MIN(Ce,0.)*grad_wb(1,j-1))/               &
     &                        (cff+Cx)

               if ((dSdt*dSdx) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
                 	inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt
               if (inflow .eq. 1) then
                   obc_sw(j,k) = sao_wb1(1,j,k) + tau *                 &
     &                           (cli_sNw(j,k) - sao_wb1(1,j,k))
               else
                   obc_sw(j,k) = obc_sw(j,k) + cli_sw(j,k)
#ifdef RADIATION_NUDGING
                   obc_sw(j,k) = obc_sw(j,k) +                  &
     &                       tau*(cli_sw(j,k)-                          &
     &                            sao_wb1(1,j,k) )

#endif /*RADIATION_NUDGING*/
               endif
               cli_sNw(j,k) = cli_sw(j,k)
              endif
            enddo
            cli_sNw(je,k) = cli_sw(je,k)
            cli_sNw(1, k) = cli_sw(1, k)
       enddo
       sao_wb1=sao_wb0
      end if

!
!---------------------------------------------------------------------------
!  Eastern edge
!---------------------------------------------------------------------------
!
      if (EAST .and. have_g_ie .and. TSTEast) then
         sao_eb0=SSS(ie2:ie,:,:)
         do k=1,ke
         	  grad_eb=0.0
            do j=je1,1,-1
               grad_eb(ie1,j) = sao_eb1(ie1,j,k) - cli_sNe(j,k) -       &
     &                          sao_eb1(ie1,j+1,k) + cli_sNe(j+1,k)
               grad_eb(ie, j) = sao_eb1(ie, j,k)  - cli_sNe(j,k) -       &
     &                          sao_eb1(ie ,j+1,k) + cli_sNe(j+1,k)
             enddo
! Attentions: the corners excluded
            do j=je1,2,-1
            	if( weto(ie,j,k) .eq. 1) then
               dSdt=sao_eb1(ie1,j,k) - cli_sNe(j,k) -                   &
     &                          sao_eb0(ie1,j,k) + cli_se(j,k)
               dSdx=sao_eb0(ie1,j,k)  - cli_se(j,k)  -                   &
     &                          sao_eb0(ie2,j,k) + cli_se(j,k)

            if ((dSdt*dSdx) .lt. 0.) dSdt=0.
            if ((dSdt*(grad_eb(ie1,j)+grad_eb(ie1,j-1))) .gt. 0.) then
               dSde = grad_eb(ie1,j)
            else
               dSde = grad_eb(ie1,j-1)
            end if
            cff = max(dSdx*dSdx+dSde*dSde,obcerr)
            Cx = dSdt*dSdx
#ifdef RADIATION_2D
            Ce = MIN(cff,MAX(dSdt*dSde,-cff))
#else
            Ce = 0.
#endif /*RADIATION_2D*/
            obc_se(j,k) = (cff*(sao_eb1(ie, j,k) - cli_sNe(j,k)) +      &
     &                     Cx *(sao_eb0(ie1,j,k) - cli_se(j,k) ) -      &
     &                     MAX(Ce,0.)*grad_eb(ie,j  )-                  &
     &                     MIN(Ce,0.)*grad_eb(ie,j-1))/                 &
     &                     (cff+Cx)

               if ((dSdt*dSdx) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
               	  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

               if (inflow .eq. 1) then
                obc_se(j,k) = sao_eb1(ie,j,k) + tau *                   &
     &                           (cli_sNe(j,k) - sao_eb1(ie,j,k))
               else
                obc_se(j,k) = obc_se(j,k) + cli_se(j,k)
#ifdef RADIATION_NUDGING
               obc_se(j,k) = obc_se(j,k) +                              &
     &                       tau*(cli_se(j,k)-                          &
     &                            sao_eb1(ie,j,k) )

#endif /*RADIATION_NUDGING*/
               endif
               cli_sNe(j,k) = cli_se(j,k)
              endif
            end do
            cli_sNe(je,k) = cli_se(je,k)
            cli_sNe(1, k) = cli_se(1, k)
         end do
         sao_eb1=sao_eb0
      endif
!----------------------------------------------------------------------------
!  Southern edge
!----------------------------------------------------------------------------
      if (SOUTH .and. have_g_je .and. TSTSouth) then
         sao_sb0=SSS(:,je2:je,:)
         do k=1,ke
         	  grad_sb=0.0
            do i=2,ie
               grad_sb(i,je1)= sao_sb1(i,je1,k) - cli_sNs(i,k) -        &
     &                         sao_sb1(i-1,je1,k) + cli_sNs(i-1,k)
               grad_sb(i,je) = sao_sb1(i,je,k) - cli_sNs(i,k) -         &
     &                         sao_sb1(i-1,je,k) + cli_sNs(i-1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
               	if( weto(i,je,k) .eq. 1) then
                  dSdt=sao_sb1(i,je1,k) - cli_sNs(i,k) -                &
     &                         sao_sb0(i,je1,k) + cli_ss(i,k)
                  dSde=sao_sb0(i,je1,k)  - cli_ss(i,k)  -                &
     &                         sao_sb0(i,je2,k) + cli_ss(i,k)

                  if ((dSdt*dSde) .lt. 0.) dSdt=0.
                  if ((dSdt*(grad_sb(i,je1)+grad_sb(i+1,je1))) .gt.0.)then
                     dSdx = grad_sb(i,je1)
                  else
                     dSdx = grad_sb(i+1,je1)
                  end if
                  cff=MAX(dSdx*dSdx+dSde*dSde,obcerr)
                  Ce=dSdt*dSde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dSdt*dSdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_ss(i,k) = (cff*(sao_sb1(i,je, k)-cli_sNs(i,k)) +  &
     &                           Ce *(sao_sb0(i,je1,k)-cli_ss(i,k) ) -  &
     &                           MAX(Cx,0.)*grad_sb(i,je)-              &
     &                           MIN(Cx,0.)*grad_sb(i+1,je))/           &
     &                           (cff+Ce)
               if ((dSdt*dSde) .lt. 0.0) then
                	inflow = 1
                  tau = tobc_in
               else
                  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

               if (inflow .eq. 1) then
                     obc_ss(i,k) = sao_sb1(i,je,k) + tau *              &
     &                             (cli_sNs(i,k) - sao_sb1(i,je,k))
                  else
                     obc_ss(i,k) = obc_ss(i,k) + cli_ss(i,k)
#ifdef RADIATION_NUDGING
                      obc_ss(i,k) = obc_ss(i,k) +                &
     &                       tau*(cli_ss(i,k)-                          &
     &                            sao_sb1(i,je,k) )

#endif /*RADIATION_NUDGING*/
                endif
                cli_sNs(i,k) = cli_ss(i,k)
               endif
               end do
               cli_sNs(ie,k) = cli_ss(ie,k)
               cli_sNs(1, k) = cli_ss(1, k)
         end do
         sao_sb1=sao_sb0
      endif

!-----------------------------------------------------------------------------
!  Northern edge
!-----------------------------------------------------------------------------

      if (NORTH .and. have_g_js .and. TSTNorth) then
         sao_nb0=SSS(:,1:3,:)
         do k=1,ke
            grad_nb=0.0
            do i=2,ie
               grad_nb(i,2) = sao_nb1(i,2,k) - cli_sNn(i,k) -           &
     &                       sao_nb1(i-1,2,k) + cli_sNn(i-1,k)
               grad_nb(i,1) = sao_nb1(i,1,k) - cli_sNn(i,k) -           &
     &                       sao_nb1(i-1,1,k) + cli_sNn(i-1,k)
            end do
! Attentions: the corners excluded
               do i=2,ie1
               	if( weto(i,1,k) .eq. 1) then
                  dSdt=sao_nb1(i,2,k) - cli_sNn(i,k) -                  &
     &                       sao_nb0(i,2,k) + cli_sn(i,k)
                  dSde=sao_nb0(i,2,k) - cli_sn(i,k) -                   &
     &                       sao_nb0(i,3,k) + cli_sn(i,k)

                  if ((dSdt*dSde) .lt. 0.) dSdt=0.
                  if ((dSdt*(grad_nb(i,2)+grad_nb(i+1,2))) .gt. 0.) then
                     dSdx = grad_nb(i,2)
                  else
                     dSdx = grad_nb(i+1,2)
                  end if
                  cff=MAX(dSdx*dSdx+dSde*dSde,obcerr)
                  Ce=dSdt*dSde
#ifdef RADIATION_2D
                  Cx = MIN(cff,MAX(dSdt*dSdx,-cff))
#else
                  Cx = 0.
#endif /*RADIATION_2D*/
                  obc_sn(i,k) = (cff*(sao_nb1(i,1,k)-cli_sNn(i,k)) +    &
     &                          Ce *(sao_nb0(i,2,k)-cli_sn(i,k)  ) -    &
     &                          MAX(Cx,0.)*grad_nb(i,1)-                &
     &                          MIN(Cx,0.)*grad_nb(i+1,1))/             &
     &                          (cff+Ce)

               if ((dSdt*dSde) .lt. 0.0) then
               	  inflow = 1
                  tau = tobc_in
               else
               	  inflow = 0
                  tau = tobc_out
               end if
               tau = tau*dt

               if (inflow .eq. 1) then
                     obc_sn(i,k) = sao_nb1(i,1,k) + tau *               &
     &                             (cli_sNn(i,k) - sao_nb1(i,1,k))
               else
                     obc_sn(i,k) = obc_sn(i,k) + cli_sn(i,k)
#ifdef RADIATION_NUDGING
                     obc_sn(i,k) = obc_sn(i,k) +                 &
     &                       tau*(cli_sn(i,k)-                          &
     &                            sao_nb1(i,1,k) )
#endif /*RADIATION_NUDGING*/
                endif
                cli_sNn(i,k) = cli_sn(i,k)
               endif
               end do
               cli_sNn(ie,k) = cli_sn(ie,k)
               cli_sNn(1, k) = cli_sn(1, k)
         end do
         sao_nb1=sao_nb0
      endif

      CALL UPDATE_ORLANSKI_TS(TTT,SSS) ! in mo_obcs.f90

      TTT0(:,:,:)=TTT(:,:,:)
      SSS0(:,:,:)=SSS(:,:,:)
#endif /*OBCS_TST_NEST*/
      END SUBROUTINE OBCS_TST_TS

! H.H --For SeaLevel Nudging
      SUBROUTINE UPGRADE_OBCS_ETA
#ifndef OBCS_UV_FLATHER
! Relaxation
      if ( West .AND. have_g_is) then
         obc_zw(:)=(obc_zwt(:)+obc_zws(:))*weto(1,:,1)
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then
         obc_ze(:)=(obc_zet(:)+obc_zes(:))*weto(ie,:,1)
      ENDIF

      ! north OB
      if ( North .AND. have_g_js) then
         obc_zn(:)=(obc_znt(:)+obc_zns(:))*weto(:,1,1)
      endif

      ! south OB
      if ( South .AND. have_g_je) then
         obc_zs(:)=(obc_zst(:)+obc_zss(:))*weto(:,je,1)
      endif
#else
!
#ifdef OBCSSHRELAX
! The relaxation OB differs from the above,
! when using tst boundary with relaxaion scheme about zeta
      if ( West .AND. have_g_is) then
         obc_zw(:)=(zeta_subT(1,:) +ZetaT(1,:))*weto(1,:,1)
        ! WRITE(IO_STDOUT,*) ' UPGRADE_OBCS_ETA West  ',obc_zw(2:10)
        ! WRITE(IO_STDOUT,*) ' UPGRADE_OBCS_ETA West  ',ZetaT(1,2:10)
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then
         obc_ze(:)=(zeta_subT(ie,:)+ZetaT(ie,:))*weto(ie,:,1)
      ENDIF

      ! north OB
      if ( North .AND. have_g_js) then
         obc_zn(:)=(zeta_subT(:,1)+ZetaT(:,1))*weto(:,1,1)
        ! WRITE(IO_STDOUT,*) ' UPGRADE_OBCS_ETA North  ',obc_zn
      endif

      ! south OB
      if ( South .AND. have_g_je) then
         obc_zs(:)=(zeta_subT(:,je)+ZetaT(:,je))*weto(:,je,1)
        ! WRITE(IO_STDOUT,*) ' UPGRADE_OBCS_ETA South  ',obc_zs
      endif
#endif
#endif
      END SUBROUTINE UPGRADE_OBCS_ETA

      SUBROUTINE UPGRADE_OBCS_UV
#ifndef OBCS_UV_FLATHER
! Pure Relaxation
      ! west OB
      if ( West .AND. have_g_is) then
       do k = 1,ke
         obc_uw(:,k)=(obc_uzw(:) + obc_usw(:,k))*amsuo(0,:,k)
         obc_vw(:,k)=(obc_vzw(:) + obc_vsw(:,k))*amsue(1,:,k)
        enddo
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then
       do k = 1,ke
         obc_ue(:,k)=(obc_uze(:) + obc_use(:,k))*amsuo(ie,:,k)
         obc_ve(:,k)=(obc_vze(:) + obc_vse(:,k))*amsue(ie,:,k)
        ENDDO
      ENDIF

      ! north OB
      if ( North .AND. have_g_js) then
       do k = 1,ke
        obc_un(:,k)=(obc_uzn(:) + obc_usn(:,k))*amsuo(:,1,k)
        obc_vn(:,k)=(obc_vzn(:) + obc_vsn(:,k))*amsue(:,0,k)
       enddo
      endif

      ! south OB
      if ( South .AND. have_g_je) then
       do k = 1,ke
        obc_us(:,k)=(obc_uzs(:) + obc_uss(:,k))*amsuo(:,je,k)
        obc_vs(:,k)=(obc_vzs(:) + obc_vss(:,k))*amsue(:,je,k)
       enddo
      endif
#else
!
#ifdef OBCUVRELAX
! The relaxation OB differs from the above,
! when using tst boundary with relaxaion scheme about UV
   if ( West .AND. have_g_is) then
       do k = 1,ke
         obc_uw(:,k)=(obc_uzw(:) + obc_usw(:,k))*amsuo(0,:,k)
         obc_vw(:,k)=(VSET(1,:) + obc_vsw(:,k))*amsue(1,:,k)
        enddo
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then
       do k = 1,ke
         obc_ue(:,k)=(obc_uze(:) + obc_use(:,k))*amsuo(ie,:,k)
         obc_ve(:,k)=(VSET(ie,:) + obc_vse(:,k))*amsue(ie,:,k)
        ENDDO
      ENDIF


      ! north OB
      if ( North .AND. have_g_js) then
       do k = 1,ke
        obc_un(:,k)=(USOT(:,1) + obc_usn(:,k))*amsuo(:,1,k)
        obc_vn(:,k)=(obc_vzn(:) + obc_vsn(:,k))*amsue(:,0,k)
       enddo
      endif

      ! south OB
      if ( South .AND. have_g_je) then
       do k = 1,ke
        obc_us(:,k)=(USOT(:,je) + obc_uss(:,k))*amsuo(:,je,k)
        obc_vs(:,k)=(obc_vzs(:) + obc_vss(:,k))*amsue(:,je,k)
       enddo
      endif
#endif
#endif

      END SUBROUTINE UPGRADE_OBCS_UV
!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!----------------------- ORCTM OBC -------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
      SUBROUTINE OCOBCS(i_step)
!----------------------------------------------------------------------------------------------------------------------------------------
!  OBCS ---- Nesting conditions
#ifdef OBCS_TST_NEST
     CALL FORC_TST_TIDES(i_step)
#endif /*OBCS_TST_NEST*/
!
!----------------------------------------------------------------------------------------------------------------------------------------
!   OBCS  ---- ZO
#ifdef OBC_ETA_FORC
     CALL FORC_OBCS_ETA(i_step) 

!   RELAXATION  BOUNDARY FORCINGS  (DEFAULT)
#ifndef OBCS_UV_FLATHER
      CALL UPGRADE_OBCS_ETA
!         Forced ETA VIA RELAXATION LAYERS
#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_ETA_EXP(ZO) 
#else
      CALL SPONGE_OBCS_ETA(ZO) 
#endif    /*OBC_SPONGE_EXP*/
!
#else ! H.H -- Flather scheme module DOBCS_UV_FLATHER
!
       ! WRITE(IO_STDOUT,*) 'Flather Boundary on/off:  on'
#ifdef OBC_TIDE
! -- Upgrade the boundary Tidal Level (2D) value.(ZetaT)
         CALL UPGRADE_FLA_BT_ZO
#endif /*OBC_TIDE*/
!
#ifdef OBC_SUBTIDE
! -- Upgrade the boundary SubTidal Level (2D) value.(zeta_subT)
         CALL UPGRADE_FLA_SUBT_ZO
#endif /*OBC_SUBTIDE*/
!
!
#ifdef OBCSSHRELAX
! --  OB Relaxation Layers settings and switch added
      CALL UPGRADE_OBCS_ETA
      CALL SPONGE_OBCS_ETA(ZO) 
#endif   /*OBCSSHRELAX*/
!
#endif   /*OBCS_UV_FLATHER*/
#endif  /*OBC_ETA_FORC*/
!
!----------------------------------------------------------------------------------------------------------------------------------------
!   OBCS  ---- UV
!   RELAXATION  BOUNDARY FORCINGS (DEFAULT)
#ifdef OBC_UV_FORC
      CALL FORC_OBCS_UV(i_step)
!
#ifndef OBCS_UV_FLATHER
      CALL UPGRADE_OBCS_UV
#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_UV_EXP(UOO,VOE)
#else
      CALL SPONGE_OBCS_UV(UOO,VOE)
#endif /*OBC_SPONGE_EXP*/
!
#endif /*OBCS_UV_FLATHER*/
!
!
!----------------------------------------------------------------------------------------------------------------------------------------
!  TIDAL AND SUBTIDAL  BOUNDARY FORCINGS ---- UV
#ifdef OBCS_UV_FLATHER
!
! --------------------------------------------------------------------------------------------------------------------------------------------------!
! ----------------------- Barotropic OBC ----------------------------------------------------------------------------------------------------!
! --------------------------------------------------------------------------------------------------------------------------------------------------!
#ifdef OBC_TIDE
! -- Upgrade the boundary Tidal velocity (2D) value.(USOT & VSET)
         CALL UPGRADE_FLA_BT_UV
         ! For FLA&ORLANSKI-OBC
#endif
!
#if defined OBC_SUBTIDE && defined ORLANSKI
! -- Upgrade the boundary Nontidal velocity (3D) value.(obcf)
        CALL UPGRADE_ORL_SUBT_UV 
        ! For FLA&ORLANSKI-OBC
#endif
!
#ifdef OBCS_TST_NEST
! -- Decomposition about the boundary Tidal and Subtidal values.(USOT, VSET, and ZetaT)
         CALL PRE_TIDES
        ! CALL SET_STL(ZO)
        ! For TST-OBC
#endif /*OBCS_TST_NEST*/
!
! --------------------------------------------------------------------------------------------------------------------------------------------------!
! ----------------------- Baroclinic OBC -----------------------------------------------------------------------------------------------------!
! --------------------------------------------------------------------------------------------------------------------------------------------------!
!
! --  OB values Closed
         CALL SET_OBCS_UV(UOO,VOE)
        ! update VOE UOO on the OB grids
!
!
#ifdef OBCUVRELAX
! --  OB Relaxation Layers settings and switch added
         CALL UPGRADE_OBCS_UV
         CALL SPONGE_OBCS_UV(UOO,VOE)
#endif
!
!
#ifdef OBCS_TST_NEST
!  3D form of Tidal and subtidal OBC basd on Liu and Gan 2016
        CALL OBCS_TST_UVK(UOO,VOE,i_step)
#else
#if defined ORLANSKI && defined OBCS_UV_FLATHER
!   3D form of Nudging or Radiation OBC based on ORLANSKI (1976)
        CALL ORLANSKI_UV(UOO,VOE,i_step)
#endif
#endif
!
#endif /*OBCS_UV_FLATHER*/
#endif /*OBC_UV_FORC*/
     
      END SUBROUTINE OCOBCS
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!----------------------- Baroclinic OBC -----------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------------------------!
      SUBROUTINE OCOBCSTS(i_step)

!----------------------------------------------------------------------------------------------------------------------------------------
!   OBCS  ---- TRANCERS(T & S)
!   RELAXATION  BOUNDARY NUDGED CLIM FORCINGS (DEFAULT)
#ifdef OBC_TS_FORC
      CALL FORC_OBCS_TS(i_step)
!
#ifndef OBCS_UV_FLATHER
#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_TS_EXP(THO,SAO) 
#else
      CALL SPONGE_OBCS_TS(THO,SAO) 
#endif /*OBC_SPONGE_EXP*/
!
#else
#ifdef OBCS_TST_NEST
!  3D form of temperature and salinity OBC basd on Liu and Gan 2016
      ! CALL OBCS_TST_TS(THO,SAO,i_step)
     
#else
      ! CALL ORLANSKI_TS(THO,SAO,i_step) ! 2-D and Nudging form
      ! CALL PRED_OBCS_TS(THO,SAO) ! 1D
#endif /*OBCS_TST_NEST*/

#ifdef OBCTSRELAX
! ! --  OB Relaxation Layers settings and switch added
      CALL SPONGE_OBCS_TS(THO,SAO) 
#endif /*OBCTSRELAX*/

#endif /*OBCS_UV_FLATHER*/

#endif /*OBC_TS_FORC*/
      END  SUBROUTINE OCOBCSTS

      SUBROUTINE OCOBCSBT(i_step)
#ifdef OBCS_TST_NEST
       DO J=1,JE
       DO I=I_start,IE
       USO(I,J)=USO(I,J)*DEUTIO(I,J)
       ENDDO
       ENDDO
!
       DO J=J_start,JE
       DO I=1,IE
       VSE(I,J)=VSE(I,J)*DEUTIE(I,J)
       ENDDO
       ENDDO
      

       CALL OBCS_TST_UVO(i_step)
     
       DO J=1,JE
       DO I=I_start,IE
       USO(I,J)=USO(I,J)*DEUTO(I,J)
       ENDDO
       ENDDO
     
       DO J=J_start,JE
       DO I=1,IE
       VSE(I,J)=VSE(I,J)*DEUTE(I,J)
       ENDDO
       ENDDO
     
#endif
      END  SUBROUTINE OCOBCSBT
     
     
     
      END MODULE MO_OBCS_NEST
