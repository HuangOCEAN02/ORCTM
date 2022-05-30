      MODULE MO_OBCS

      USE mo_commo1
      USE mo_param1
      USE mo_mpi
      USE mo_parallel
      USE mo_units
! Main Switch for Open Boundary 
      LOGICAL:: West = .FALSE., East = .FALSE.
      LOGICAL:: North = .FALSE., South = .FALSE.
! Radiation Condition Switch for BC flow or Total flow      
      LOGICAL:: OrlanskiWest = .FALSE., OrlanskiEast = .FALSE.
      LOGICAL:: OrlanskiNorth = .FALSE., OrlanskiSouth = .FALSE.
! Sponge Forcing Condition Switch for BT flow or Elevation
      LOGICAL:: SpongeWest = .FALSE., SpongeEast = .FALSE.
      LOGICAL:: SpongeNorth = .FALSE., SpongeSouth = .FALSE. 
! whether tidal wave applies to the 2D-Sponge Forcing Condition(BT flow)
      LOGICAL::add_tide = .FALSE.     
! Partial tide numbers             
      INTEGER,PARAMETER::tidalnumber = 2     
! TS-Read Time Interval（units: [s]）                  
      REAL,PARAMETER::cli_step=86400          
! Orlanski radiation condition parameters         
      REAL,PARAMETER::obcerr=1.E-10,ab1=1.51,ab2=-0.51, &
                    time_scale=43200,obrel_tem=5.E-8,obrel_sal=5.E-8
! Sponge layer grid numbers on each side                
      INTEGER::spongew,spongee,spongen,sponges 
! H.H-- (1/x) Sponge layer time parameters(TS,UV & ETA,units: [s])      
      REAL::Urelaxinner,Urelaxbound,Vrelaxinner,Vrelaxbound  
! H.H-- (Exp) Sponge layer time parameters(TS,UV & ETA,units: [s])        
      REAL::TAO                      
! Sponge layer time parameters                     
      REAL:: rel_tsi,rel_tsb,rel_uvi,rel_uvb,rel_etai,rel_etab 
! Sponge layer thickness
      REAL,POINTER::Lspp_W(:),Lspp_E(:),Lspp_N(:),Lspp_S(:)    
      REAL,POINTER::Lspu_W(:),Lspu_E(:),Lspu_N(:),Lspu_S(:)   
      REAL,POINTER::Lspv_W(:),Lspv_E(:),Lspv_N(:),Lspv_S(:)
! tidal period                       
      REAL,POINTER :: tidalPeriod(:)
! tidal harmonic constants     
      REAL,POINTER :: OBWZam(:,:),OBEZam(:,:),OBNZam(:,:),OBSZam(:,:), &
     &         OBWZph(:,:),OBEZph(:,:),OBNZph(:,:),OBSZph(:,:),        &
     &         OBWUam(:,:),OBEUam(:,:),OBNUam(:,:),OBSUam(:,:),        &
     &         OBWUph(:,:),OBEUph(:,:),OBNUph(:,:),OBSUph(:,:),        &
     &         OBWVam(:,:),OBEVam(:,:),OBNVam(:,:),OBSVam(:,:),        &
     &         OBWVph(:,:),OBEVph(:,:),OBNVph(:,:),OBSVph(:,:)         
     
     
      REAL, POINTER :: ZFO(:,:),ZFO_W(:,:),ZFO_E(:,:),ZFO_N(:,:),ZFO_S(:,:), &
     &         UFO(:,:),UFO_W(:,:),UFO_E(:,:),UFO_N(:,:),UFO_S(:,:),    &
     &         VFE(:,:),VFE_W(:,:),VFE_E(:,:),VFE_N(:,:),VFE_S(:,:)
     
      REAL, POINTER :: UFKO(:,:,:),UFKO_W(:,:,:),UFKO_E(:,:,:),UFKO_N(:,:,:),UFKO_S(:,:,:), &
     &         VFKE(:,:,:),VFKE_W(:,:,:),VFKE_E(:,:,:),VFKE_N(:,:,:),VFKE_S(:,:,:),    &
     &         TFKO(:,:,:),TFKO_W(:,:,:),TFKO_E(:,:,:),TFKO_N(:,:,:),TFKO_S(:,:,:),    &
     &         SFKO(:,:,:),SFKO_W(:,:,:),SFKO_E(:,:,:),SFKO_N(:,:,:),SFKO_S(:,:,:)
     
      REAL, POINTER :: uoo_wb0(:,:,:),uoo_wb1(:,:,:),uoo_wb2(:,:,:),    &
     &                 voe_wb0(:,:,:),voe_wb1(:,:,:),voe_wb2(:,:,:),    &
     &                 uoo_eb0(:,:,:),uoo_eb1(:,:,:),uoo_eb2(:,:,:),    &
     &                 voe_eb0(:,:,:),voe_eb1(:,:,:),voe_eb2(:,:,:),    &
     &                 uoo_nb0(:,:,:),uoo_nb1(:,:,:),uoo_nb2(:,:,:),    &
     &                 voe_nb0(:,:,:),voe_nb1(:,:,:),voe_nb2(:,:,:),    &
     &                 uoo_sb0(:,:,:),uoo_sb1(:,:,:),uoo_sb2(:,:,:),    &
     &                 voe_sb0(:,:,:),voe_sb1(:,:,:),voe_sb2(:,:,:)

      REAL, POINTER :: tho_wb0(:,:,:),tho_wb1(:,:,:),tho_wb2(:,:,:),    &
     &                 sao_wb0(:,:,:),sao_wb1(:,:,:),sao_wb2(:,:,:),    &
     &                 tho_eb0(:,:,:),tho_eb1(:,:,:),tho_eb2(:,:,:),    &
     &                 sao_eb0(:,:,:),sao_eb1(:,:,:),sao_eb2(:,:,:),    &
     &                 tho_nb0(:,:,:),tho_nb1(:,:,:),tho_nb2(:,:,:),    &
     &                 sao_nb0(:,:,:),sao_nb1(:,:,:),sao_nb2(:,:,:),    &
     &                 tho_sb0(:,:,:),tho_sb1(:,:,:),tho_sb2(:,:,:),    &
     &                 sao_sb0(:,:,:),sao_sb1(:,:,:),sao_sb2(:,:,:)

      REAL, POINTER :: obc_zw(:),obc_uw(:,:),obc_vw(:,:),obc_uzw(:),obc_vzw(:), &
     &         obc_ze(:),obc_ue(:,:),obc_ve(:,:),obc_uze(:),obc_vze(:), &
     &         obc_zn(:),obc_un(:,:),obc_vn(:,:),obc_uzn(:),obc_vzn(:), &
     &         obc_zs(:),obc_us(:,:),obc_vs(:,:),obc_uzs(:),obc_vzs(:)
     
     
      REAL, POINTER :: obc_tw(:,:),obc_sw(:,:),cli_tw(:,:),cli_sw(:,:), &
     &                 obc_te(:,:),obc_se(:,:),cli_te(:,:),cli_se(:,:), &
     &                 obc_tn(:,:),obc_sn(:,:),cli_tn(:,:),cli_sn(:,:), &
     &                 obc_ts(:,:),obc_ss(:,:),cli_ts(:,:),cli_ss(:,:)

      REAL, POINTER :: cl_uw(:,:),cl_vw(:,:),cl_tw(:,:),cl_sw(:,:),     &
     &                 cl_ue(:,:),cl_ve(:,:),cl_te(:,:),cl_se(:,:),     &
     &                 cl_un(:,:),cl_vn(:,:),cl_tn(:,:),cl_sn(:,:),     &
     &                 cl_us(:,:),cl_vs(:,:),cl_ts(:,:),cl_ss(:,:)

      INTEGER,PARAMETER::                                               &
     & IO_IN_OBCZW=911,IO_IN_OBCZE=921,IO_IN_OBCZN=931,IO_IN_OBCZS=941, &
     & IO_IN_OBCTW=912,IO_IN_OBCTE=922,IO_IN_OBCTN=932,IO_IN_OBCTS=942, &
     & IO_IN_OBCSW=915,IO_IN_OBCSE=925,IO_IN_OBCSN=935,IO_IN_OBCSS=945, &
     & IO_IN_OBCUW=913,IO_IN_OBCUE=923,IO_IN_OBCUN=933,IO_IN_OBCUS=943, &
     & IO_IN_OBCVW=914,IO_IN_OBCVE=924,IO_IN_OBCVN=934,IO_IN_OBCVS=944, &
     & IO_IN_OBCZGW=811,IO_IN_OBCZGE=821,IO_IN_OBCZGN=831,IO_IN_OBCZGS=841, &
     & IO_IN_OBCUGW=813,IO_IN_OBCUGE=823,IO_IN_OBCUGN=833,IO_IN_OBCUGS=843, &
     & IO_IN_OBCVGW=814,IO_IN_OBCVGE=824,IO_IN_OBCVGN=834,IO_IN_OBCVGS=844
     
     
     


      CONTAINS


! allocate memory for OBCS module
      SUBROUTINE alloc_mem_obcs

     ! components   = M2       S2     N2       K2       K1      O1       P1       Q1      
     ! periods (hr) = 12.4206  12     12.6583  11.9672  23.9345 25.8193  24.0659  26.8684  
     allocate(tidalPeriod(1:tidalnumber))         
     allocate( Lspp_W(JE),Lspp_E(JE),Lspp_N(IE),Lspp_S(IE),               &
     &         Lspu_W(JE),Lspu_E(JE),Lspu_N(0:IE),Lspu_S(0:IE),           &
     &         Lspv_W(0:JE),Lspv_E(0:JE),Lspv_N(IE),Lspv_S(IE))
     allocate( OBWZam(1:je_g,1:tidalnumber),OBEZam(1:je_g,1:tidalnumber), &
     &         OBNZam(1:ie_g,1:tidalnumber),OBSZam(1:ie_g,1:tidalnumber), &
     &         OBWZph(1:je_g,1:tidalnumber),OBEZph(1:je_g,1:tidalnumber), &
     &         OBNZph(1:ie_g,1:tidalnumber),OBSZph(1:ie_g,1:tidalnumber), &
     
     &         OBWUam(1:je_g,1:tidalnumber),OBEUam(1:je_g,1:tidalnumber), &
     &         OBNUam(0:ie_g,1:tidalnumber),OBSUam(0:ie_g,1:tidalnumber), &
     &         OBWUph(1:je_g,1:tidalnumber),OBEUph(1:je_g,1:tidalnumber), &
     &         OBNUph(0:ie_g,1:tidalnumber),OBSUph(0:ie_g,1:tidalnumber), &
     
     &         OBWVam(0:je_g,1:tidalnumber),OBEVam(0:je_g,1:tidalnumber), &
     &         OBNVam(1:ie_g,1:tidalnumber),OBSVam(1:ie_g,1:tidalnumber), &
     &         OBWVph(0:je_g,1:tidalnumber),OBEVph(0:je_g,1:tidalnumber), &
     &         OBNVph(1:ie_g,1:tidalnumber),OBSVph(1:ie_g,1:tidalnumber))


      allocate( ZFO(ie,je),  UFO(I_start:ie,je),  VFE(ie,J_start:je),   &
     &          ZFO_W(ie,je),UFO_W(I_start:ie,je),VFE_W(ie,J_start:je), &
     &          ZFO_E(ie,je),UFO_E(I_start:ie,je),VFE_E(ie,J_start:je), &
     &          ZFO_N(ie,je),UFO_N(I_start:ie,je),VFE_N(ie,J_start:je), &
     &          ZFO_S(ie,je),UFO_S(I_start:ie,je),VFE_S(ie,J_start:je))

      allocate( UFKO(I_start:ie,je,ke),VFKE(ie,J_start:je,ke),   &
     &          UFKO_W(I_start:ie,je,ke),VFKE_W(ie,J_start:je,ke), &
     &          UFKO_E(I_start:ie,je,ke),VFKE_E(ie,J_start:je,ke), &
     &          UFKO_N(I_start:ie,je,ke),VFKE_N(ie,J_start:je,ke), &
     &          UFKO_S(I_start:ie,je,ke),VFKE_S(ie,J_start:je,ke))
     
      allocate( TFKO(ie,je,ke),SFKO(ie,je,ke),  &
     &          TFKO_W(ie,je,ke),SFKO_W(ie,je,ke), &
     &          TFKO_E(ie,je,ke),SFKO_E(ie,je,ke), &
     &          TFKO_N(ie,je,ke),SFKO_N(ie,je,ke), &
     &          TFKO_S(ie,je,ke),SFKO_S(ie,je,ke)) 

      allocate( obc_zw(je) )
      allocate( uoo_wb0(0:2,je,ke), voe_wb0(1:3,J_start:je,ke),         &
     &          uoo_wb1(0:2,je,ke), voe_wb1(1:3,J_start:je,ke),         &
     &          uoo_wb2(0:2,je,ke), voe_wb2(1:3,J_start:je,ke),         &
     &               obc_uw(je,ke),      obc_vw(J_start:je,ke),         &
     &              obc_uzw(je   ),     obc_vzw(J_start:je   ),         &
     &                cl_uw(je,ke),       cl_vw(J_start:je,ke) )         
     
     
      allocate( tho_wb0(1:3,je,ke), sao_wb0(1:3,je,ke),                 &
     &          tho_wb1(1:3,je,ke), sao_wb1(1:3,je,ke),                 &
     &          tho_wb2(1:3,je,ke), sao_wb2(1:3,je,ke),                 &
     &               obc_tw(je,ke),      obc_sw(je,ke),                 &
     &               cli_tw(je,ke),      cli_sw(je,ke),                 &
     &                cl_tw(je,ke),       cl_sw(je,ke) )

      allocate( obc_ze(je) )
      allocate( uoo_eb0(ie2:ie,je,ke), voe_eb0(ie2:ie,J_start:je,ke),   &
     &          uoo_eb1(ie2:ie,je,ke), voe_eb1(ie2:ie,J_start:je,ke),   &
     &          uoo_eb2(ie2:ie,je,ke), voe_eb2(ie2:ie,J_start:je,ke),   &
     &                  obc_ue(je,ke),         obc_ve(J_start:je,ke),   &
     &                 obc_uze(je   ),        obc_vze(J_start:je   ),   &
     &                   cl_ue(je,ke),          cl_ve(J_start:je,ke) )
     
          
      allocate( tho_eb0(ie2:ie,je,ke), sao_eb0(ie2:ie,je,ke),           &
     &          tho_eb1(ie2:ie,je,ke), sao_eb1(ie2:ie,je,ke),           &
     &          tho_eb2(ie2:ie,je,ke), sao_eb2(ie2:ie,je,ke),           &
     &                  obc_te(je,ke),         obc_se(je,ke),           &
     &                  cli_te(je,ke),         cli_se(je,ke),           &
     &                   cl_te(je,ke),          cl_se(je,ke) )

      allocate( obc_zn(ie) )
      allocate( uoo_nb0(I_start:ie,1:3,ke), voe_nb0(ie,0:2,ke),         &
     &          uoo_nb1(I_start:ie,1:3,ke), voe_nb1(ie,0:2,ke),         &
     &          uoo_nb2(I_start:ie,1:3,ke), voe_nb2(ie,0:2,ke),         &
     &           obc_un(I_start:ie,    ke),  obc_vn(ie,    ke),         &
     &          obc_uzn(I_start:ie       ), obc_vzn(ie       ),         &
     &            cl_un(I_start:ie,    ke),   cl_vn(ie,    ke) )
     
      allocate( tho_nb0(ie,1:3,ke), sao_nb0(ie,1:3,ke),                 &
     &          tho_nb1(ie,1:3,ke), sao_nb1(ie,1:3,ke),                 &
     &          tho_nb2(ie,1:3,ke), sao_nb2(ie,1:3,ke),                 &
     &           obc_tn(ie,    ke),  obc_sn(ie,    ke),                 &
     &           cli_tn(ie,    ke),  cli_sn(ie,    ke),                 &
     &            cl_tn(ie,    ke),   cl_sn(ie,    ke) )

      allocate( obc_zs(ie) )
      allocate( uoo_sb0(I_start:ie,je2:je,ke), voe_sb0(ie,je2:je,ke),   &
     &          uoo_sb1(I_start:ie,je2:je,ke), voe_sb1(ie,je2:je,ke),   &
     &          uoo_sb2(I_start:ie,je2:je,ke), voe_sb2(ie,je2:je,ke),   &
     &           obc_us(I_start:ie,       ke),  obc_vs(ie,       ke),   &
     &          obc_uzs(I_start:ie          ), obc_vzs(ie          ),   &
     &            cl_us(I_start:ie,       ke),   cl_vs(ie,       ke) )          

      allocate( tho_sb0(ie,je2:je,ke), sao_sb0(ie,je2:je,ke),           &
     &          tho_sb1(ie,je2:je,ke), sao_sb1(ie,je2:je,ke),           &
     &          tho_sb2(ie,je2:je,ke), sao_sb2(ie,je2:je,ke),           &
     &           obc_ts(ie,       ke),  obc_ss(ie,       ke),           &
     &           cli_ts(ie,       ke),  cli_ss(ie,       ke),           &
     &            cl_ts(ie,       ke),   cl_ss(ie,       ke) )


      tidalPeriod(:)=0.0
      Lspp_W(:)=0.0;Lspu_W(:)=0.0;Lspv_W(:)=0.0      
      Lspp_E(:)=0.0;Lspu_E(:)=0.0;Lspv_E(:)=0.0           
      Lspp_N(:)=0.0;Lspu_N(:)=0.0;Lspv_N(:)=0.0  
      Lspp_S(:)=0.0;Lspu_S(:)=0.0;Lspv_S(:)=0.0
     
      rel_tsi=0.0;rel_tsb=0.0;rel_uvi=0.0
      rel_uvb=0.0;rel_etai=0.0;rel_etab=0.0 

      OBWZam=0.0;OBEZam=0.0;OBNZam=0.0;OBSZam=0.0
      OBWZph=0.0;OBEZph=0.0;OBNZph=0.0;OBSZph=0.0
      OBWUam=0.0;OBEUam=0.0;OBNUam=0.0;OBSUam=0.0
      OBWUph=0.0;OBEUph=0.0;OBNUph=0.0;OBSUph=0.0
      OBWVam=0.0;OBEVam=0.0;OBNVam=0.0;OBSVam=0.0
      OBWVph=0.0;OBEVph=0.0;OBNVph=0.0;OBSVph=0.0

      obc_zw=0.0
      uoo_wb0=0.0;uoo_wb1=0.0;uoo_wb2=0.0;obc_uw=0.0;cl_uw=0.0
      voe_wb0=0.0;voe_wb1=0.0;voe_wb2=0.0;obc_vw=0.0;cl_vw=0.0
      obc_uzw=0.0;obc_vzw=0.0;  
      tho_wb0=0.0;tho_wb1=0.0;tho_wb2=0.0;obc_tw=0.0;cl_tw=0.0
      sao_wb0=0.0;sao_wb1=0.0;sao_wb2=0.0;obc_sw=0.0;cl_sw=0.0

      obc_ze=0.0
      uoo_eb0=0.0;uoo_eb1=0.0;uoo_eb2=0.0;obc_ue=0.0;cl_ue=0.0
      voe_eb0=0.0;voe_eb1=0.0;voe_eb2=0.0;obc_ve=0.0;cl_ve=0.0      
      obc_uze=0.0;obc_vze=0.0;           
      tho_eb0=0.0;tho_eb1=0.0;tho_eb2=0.0;obc_te=0.0;cl_te=0.0
      sao_eb0=0.0;sao_eb1=0.0;sao_eb2=0.0;obc_se=0.0;cl_se=0.0

      obc_zn=0.0
      uoo_nb0=0.0;uoo_nb1=0.0;uoo_nb2=0.0;obc_un=0.0;cl_un=0.0
      voe_nb0=0.0;voe_nb1=0.0;voe_nb2=0.0;obc_vn=0.0;cl_vn=0.0
      obc_uzn=0.0;obc_vzn=0.0;    
      tho_nb0=0.0;tho_nb1=0.0;tho_nb2=0.0;obc_tn=0.0;cl_tn=0.0
      sao_nb0=0.0;sao_nb1=0.0;sao_nb2=0.0;obc_sn=0.0;cl_sn=0.0

      obc_zs=0.0
      uoo_sb0=0.0;uoo_sb1=0.0;uoo_sb2=0.0;obc_us=0.0;cl_us=0.0
      voe_sb0=0.0;voe_sb1=0.0;voe_sb2=0.0;obc_vs=0.0;cl_vs=0.0
      obc_uzs=0.0;obc_vzs=0.0;
      tho_sb0=0.0;tho_sb1=0.0;tho_sb2=0.0;obc_ts=0.0;cl_ts=0.0
      sao_sb0=0.0;sao_sb1=0.0;sao_sb2=0.0;obc_ss=0.0;cl_ss=0.0

      END SUBROUTINE alloc_mem_obcs


! initialize OBC forcing as model fields at time zero
      SUBROUTINE INI_OBCS
      
! OBC - Sea Surface height  
      IF(p_pe==p_io) THEN

#ifdef OBC_ETA_FORC
#ifdef OBC_TIDE
       IF( West ) THEN
        OPEN(IO_IN_OBCZW,FILE='OB_West_z_Amp.bin',  & 
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCZGW,FILE='OB_West_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
     
        READ(IO_IN_OBCZW)   ((OBWZam(i,j),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCZGW)  ((OBWZph(i,j),i=1,je_g),j=1,tidalnumber)
       ENDIF
       
       IF( East ) THEN
        OPEN(IO_IN_OBCZE,FILE='OB_East_z_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCZGE,FILE='OB_East_z_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        READ(IO_IN_OBCZE)   ((OBEZam(i,j),i=1,je_g),j=1,tidalnumber)     
        READ(IO_IN_OBCZGE)  ((OBEZph(i,j),i=1,je_g),j=1,tidalnumber)
       ENDIF
       
       IF( North ) THEN
        OPEN(IO_IN_OBCZN,FILE='OB_North_z_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCZGN,FILE='OB_North_z_Gph.bin',&
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        READ(IO_IN_OBCZN)   ((OBNZam(i,j),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_OBCZGN)  ((OBNZph(i,j),i=1,ie_g),j=1,tidalnumber)
       ENDIF

       IF( South ) THEN       
        OPEN(IO_IN_OBCZS,FILE='OB_South_z_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCZGS,FILE='OB_South_z_Gph.bin',&
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY') 
        READ(IO_IN_OBCZS)   ((OBSZam(i,j),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_OBCZGS)  ((OBSZph(i,j),i=1,ie_g),j=1,tidalnumber)         
       ENDIF
#endif      
#endif
      ENDIF
      IF( West ) THEN
        CALL p_bcast(OBWZam,p_io)
        CALL p_bcast(OBWZph,p_io)
      ENDIF
      IF( East ) THEN
        CALL p_bcast(OBEZam,p_io)
        CALL p_bcast(OBEZph,p_io)
      ENDIF
      IF( North ) THEN
        CALL p_bcast(OBNZam,p_io)
        CALL p_bcast(OBNZph,p_io)
      ENDIF
      IF( South ) THEN
        CALL p_bcast(OBSZam,p_io) 
        CALL p_bcast(OBSZph,p_io)  
      ENDIF

     
! OBC - Zonal&Meridional velocity  
      IF(p_pe==p_io) THEN
#ifdef OBC_UV_FORC
#ifdef OBC_TIDE
       IF( West ) THEN
        OPEN(IO_IN_OBCUW,FILE='OB_West_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCUGW,FILE='OB_West_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')     
        OPEN(IO_IN_OBCVW,FILE='OB_West_v_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCVGW,FILE='OB_West_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        READ(IO_IN_OBCUW)   ((OBWUam(i,j),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCUGW)  ((OBWUph(i,j),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCVW)   ((OBWVam(i,j),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCVGW)  ((OBWVph(i,j),i=0,je_g),j=1,tidalnumber)
       ENDIF
       
       IF( East ) THEN
        OPEN(IO_IN_OBCUE,FILE='OB_East_u_Amp.bin',  &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCUGE,FILE='OB_East_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')     
        OPEN(IO_IN_OBCVE,FILE='OB_East_v_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCVGE,FILE='OB_East_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        READ(IO_IN_OBCUE)   ((OBEUam(i,j),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCUGE)  ((OBEUph(i,j),i=1,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCVE)   ((OBEVam(i,j),i=0,je_g),j=1,tidalnumber)
        READ(IO_IN_OBCVGE)  ((OBEVph(i,j),i=0,je_g),j=1,tidalnumber)          
       ENDIF   
        
       IF( North ) THEN   
        OPEN(IO_IN_OBCUN,FILE='OB_North_u_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCUGN,FILE='OB_North_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCVN,FILE='OB_North_v_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCVGN,FILE='OB_North_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        READ(IO_IN_OBCUN)   ((OBNUam(i,j),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_OBCUGN)  ((OBNUph(i,j),i=0,ie_g),j=1,tidalnumber)
        READ(IO_IN_OBCVN)   ((OBNVam(i,j),i=1,ie_g),j=1,tidalnumber)
        READ(IO_IN_OBCVGN)  ((OBNVph(i,j),i=1,ie_g),j=1,tidalnumber)
       ENDIF
       
       IF( South ) THEN                  
        OPEN(IO_IN_OBCUS,FILE='OB_South_u_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCUGS,FILE='OB_South_u_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')    
        OPEN(IO_IN_OBCVS,FILE='OB_South_v_Amp.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')    
        OPEN(IO_IN_OBCVGS,FILE='OB_South_v_Gph.bin', &
     &        ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')    
        READ(IO_IN_OBCUS)   ((OBSUam(i,j),i=0,ie_g),j=1,tidalnumber)             
        READ(IO_IN_OBCUGS)  ((OBSUph(i,j),i=0,ie_g),j=1,tidalnumber)      
        READ(IO_IN_OBCVS)   ((OBSVam(i,j),i=1,ie_g),j=1,tidalnumber)     
        READ(IO_IN_OBCVGS)  ((OBSVph(i,j),i=1,ie_g),j=1,tidalnumber)     
       ENDIF    
#endif     
#endif
      ENDIF
      
      IF( West ) THEN
       CALL p_bcast(OBWUam,p_io) 
       CALL p_bcast(OBWUph,p_io)
       CALL p_bcast(OBWVam,p_io)
       CALL p_bcast(OBWVph,p_io)
      ENDIF 
      
      IF( East ) THEN                  
       CALL p_bcast(OBEUam,p_io)
       CALL p_bcast(OBEUph,p_io)
       CALL p_bcast(OBEVam,p_io)
       CALL p_bcast(OBEVph,p_io)
      ENDIF
      
      IF( North ) THEN             
       CALL p_bcast(OBNUam,p_io)
       CALL p_bcast(OBNUph,p_io)
       CALL p_bcast(OBNVam,p_io) 
       CALL p_bcast(OBNVph,p_io)
      ENDIF
       
      IF( South ) THEN                 
       CALL p_bcast(OBSUam,p_io)      
       CALL p_bcast(OBSUph,p_io)       
       CALL p_bcast(OBSVam,p_io) 
       CALL p_bcast(OBSVph,p_io)
      ENDIF
! OBC-Temperature & Salinity
      IF(p_pe==p_io) THEN
#ifdef OBC_TS_FORC

       IF( West ) THEN
        OPEN(IO_IN_OBCTW,FILE='OB_West_T.bin',                          &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCSW,FILE='OB_West_S.bin',                          &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF
       
       IF( East ) THEN        
        OPEN(IO_IN_OBCTE,FILE='OB_East_T.bin',                          &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCSE,FILE='OB_East_S.bin',                          &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF
       
      IF( North ) THEN       	           
        OPEN(IO_IN_OBCTN,FILE='OB_North_T.bin',                         &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCSN,FILE='OB_North_S.bin',                         &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
       ENDIF
       
      IF( South ) THEN                 
        OPEN(IO_IN_OBCTS,FILE='OB_South_T.bin',                         &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
        OPEN(IO_IN_OBCSS,FILE='OB_South_S.bin',                         &
     &                 ACCESS='SEQUENTIAL',ACTION='READ',FORM='BINARY')
      ENDIF
#endif
      ENDIF

      ! west OB
      if ( West ) then
      if (have_g_is) then
        obc_zw(:)  =zo(1,:)
        obc_uw(:,:)=uoo(0,:,:)
        obc_vw(:,:)=voe(1,:,:)
        obc_tw(:,:)=tho(1,:,:)
        obc_sw(:,:)=sao(1,:,:)
      endif
      endif

      ! east OB
      IF ( East ) THEN      
      if (have_g_ie) then
        obc_ze(:)  =zo(ie,:)
        obc_ue(:,:)=uoo(ie,:,:)
        obc_ve(:,:)=voe(ie,:,:)
        obc_te(:,:)=tho(ie,:,:)
        obc_se(:,:)=sao(ie,:,:)
      endif
      endif
      
      ! north OB
      IF ( North ) THEN      
      if (have_g_js) then
        obc_zn(:)  =zo(:,1)
        obc_un(:,:)=uoo(:,1,:)
        obc_vn(:,:)=voe(:,0,:)
        obc_tn(:,:)=tho(:,1,:)
        obc_sn(:,:)=sao(:,1,:)
      endif
      endif
      
      ! south OB
      IF ( South ) THEN        
      if (have_g_je) then
        obc_zs(:)  =zo(:,je)
        obc_us(:,:)=uoo(:,je,:)
        obc_vs(:,:)=voe(:,je,:)
        obc_ts(:,:)=tho(:,je,:)
        obc_ss(:,:)=sao(:,je,:)
      endif
      endif
#ifndef OBC_SPONGE_EXP
      rel_uvb  = Urelaxbound
      rel_uvi  = Urelaxinner
          
      rel_tsb  = rel_uvb
      rel_tsi  = rel_uvi
                 
      rel_etab = rel_tsb
      rel_etai = rel_tsi      
#else
      Lspp_W(:) = SUM(DLXP(1:spongew,:),DIM=1)
      Lspu_W(:) = SUM(DLXU(I_start:spongew-1,:),DIM=1)
      Lspv_W(:) = SUM(DLXV(1:spongew,:),DIM=1)
        
      Lspp_E(:) = SUM(DLXP(ie-spongee+1:ie,:),DIM=1)
      Lspu_E(:) = SUM(DLXU(ie-spongee+1:ie,:),DIM=1)
      Lspv_E(:) = SUM(DLXV(ie-spongee+1:ie,:),DIM=1)
           
      Lspp_N(:) = SUM(DLYP(:,1:spongen),DIM=2)           
      Lspu_N(:) = SUM(DLYU(:,1:spongen),DIM=2) 
      Lspv_N(:) = SUM(DLYV(:,J_start:spongen-1),DIM=2)  
           
      Lspp_S(:) = SUM(DLYP(:,je-sponges+1:je),DIM=2)            
      Lspu_S(:) = SUM(DLYU(:,je-sponges+1:je),DIM=2) 
      Lspv_S(:) = SUM(DLYV(:,je-sponges+1:je),DIM=2)      
      
      ! west OB
      if ( West ) then
       if (have_g_is) then
        write(IO_STDOUT,*) 'West Exp Sponge Layer Thickness'
        write(IO_STDOUT,*) 'Lspp_W ', Lspp_W(:)
        write(IO_STDOUT,*) 'Lspu_W ', Lspu_W(:)
        write(IO_STDOUT,*) 'Lspv_W ', Lspv_W(:) 
       endif
      endif  
     
      ! East OB
      if ( East ) then
       if (have_g_ie) then
        write(IO_STDOUT,*) 'East Exp Sponge Layer Thickness'
        write(IO_STDOUT,*) 'Lspp_E ', Lspp_E(:)
        write(IO_STDOUT,*) 'Lspu_E ', Lspu_E(:)
        write(IO_STDOUT,*) 'Lspv_E ', Lspv_E(:) 
       endif
      endif  

      ! North OB
      if ( North ) then
       if (have_g_js) then
        write(IO_STDOUT,*) 'North Exp Sponge Layer Thickness'
        write(IO_STDOUT,*) 'Lspp_N ', Lspp_N(:)
        write(IO_STDOUT,*) 'Lspu_N ', Lspu_N(:)
        write(IO_STDOUT,*) 'Lspv_N ', Lspv_N(:) 
       endif
      endif  
     
      ! South OB
      if ( South ) then
       if (have_g_je) then
        write(IO_STDOUT,*) 'South Exp Sponge Layer Thickness'
        write(IO_STDOUT,*) 'Lspp_S ', Lspp_S(:)
        write(IO_STDOUT,*) 'Lspu_S ', Lspu_S(:)
        write(IO_STDOUT,*) 'Lspv_S ', Lspv_S(:) 
       endif
      endif 
#endif 
                                              
      END SUBROUTINE INI_OBCS


! OBC eta forcing
      SUBROUTINE FORC_OBCS_ETA(i_step)

      REAL :: mytime,obc_zwt(je),obc_zet(je),obc_znt(ie),obc_zst(ie)
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
        obc_zw(:) = obc_zwt(:)
       endif
      endif
      !if ( West .AND. have_g_is ) then
        ! test for NaN in occlit.f90
      !   write(IO_STDOUT,*) 'FORC_OBCS_ETA'
      !   write(IO_STDOUT,*) 'LDTRUN', i_step
      !   write(IO_STDOUT,*) 'obc_zw', obc_zw(:)
      !endif
      
      ! east OB
      IF ( East ) THEN       
       if (have_g_ie) then
        obc_zet(:) = 0.
        DO td=1,tidalnumber
            obc_zet(:) = obc_zet(:) + OBEZam(p_joff+1:p_joff+je,td)*weto(ie,:,1)* &
             &      COS(2.0 * PI * (myTime-OBEZph(p_joff+1:p_joff+je,td))/tidalPeriod(td)) 
        ENDDO
        obc_ze(:) = obc_zet(:)
       endif
      endif
      
      !if ( East .AND. have_g_ie ) then        
        ! test for NaN in occlit.f90
      !   write(IO_STDOUT,*) 'FORC_OBCS_ETA'
      !   write(IO_STDOUT,*) 'LDTRUN', i_step
      !   write(IO_STDOUT,*) 'obc_ze', obc_ze(:)       
      !endif
      
            
      ! north OB
      IF ( North ) THEN       
       if (have_g_js) then
        obc_znt(:) = 0.
        DO td=1,tidalnumber
           obc_znt(:) = obc_znt(:) + OBNZam(p_ioff+1:p_ioff+ie,td)*weto(:,1,1)*  &
             &      COS(2.0 * PI * (myTime-OBNZph(p_ioff+1:p_ioff+ie,td))/tidalPeriod(td))
        ENDDO
        obc_zn(:) = obc_znt(:)
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
        obc_zs(:) = obc_zst(:)
       endif
      endif

#endif     
     

      END SUBROUTINE FORC_OBCS_ETA


!  Open Boundary Radiation Condition for baroclinc flow   
!  or total flow. Calculate future boundary data at the   
! open boundaries at time = futureTime by applying        
!  Orlanski radiation conditions.
! H.H -- Orlanski radiation conditions will upgrade the
! boundary values(obc_uw,obc_ue,....)
  
      SUBROUTINE PRED_OBCS_UV(UUU0,VVV0)
      REAL UUU0(1:IE-I_start+1,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,1:JE-J_start+1,KE),VVV(IE,J_start:JE,KE)

      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)
      
!  Added code to allow filtering of phase speed following
!  Blumberg and Kantha. There is now a separate array
!  cl_**, where **=Variable(U,V,T,S,W)Boundary(E,W,N,S) for
!  the dimensional phase speed. These arrays are initialized to
!  zero. CVEL_** is filtered according to
!  cl_** = f1*cl_**(new)+(1-f1)*cl_**(old).
!  f1=1.0 turns off filtering.
      f1 = dt/time_scale
      f2 = 1.0-f1
!  Changed code to average phase speed. A new variable
!  'time_scale' was created. This variable must now be
!  specified. Then, f1=dT/time_scale.Since the goal is 
!  to smooth out the 'singularities' in thediagnosed phase 
!  speed, time_scale could be picked as the duration of the
!  singular period in the unfiltered case. Thus,for a plane wave
!  time_scale might be the time take for thewave to travel 
!  a distance DX, where DX is the width of the region near which 
!  d(phi)/dx is small.  

      ! west OB
      if ( West .AND. have_g_is) then
       if (OrlanskiWest) then            
          uoo_wb0=UUU(0:2,:,:)
          do k=1,ke
            do j=1,je
              d01 = ab1*(uoo_wb0(1,J,K)-uoo_wb0(0,J,K)) +               &
     &              ab2*(uoo_wb1(1,J,K)-uoo_wb1(0,J,K))
              d12 = ab1*(uoo_wb1(2,J,K)-uoo_wb1(1,J,K)) +               &
     &              ab2*(uoo_wb2(2,J,K)-uoo_wb2(1,J,K))

              IF (abs(d12) .lt. obcerr) THEN
                 CL=0.
              ELSE
                 CL=(uoo_wb0(1,J,K)-uoo_wb1(1,J,K)) / d12
              ENDIF
              CL=min(max(CL,0.0),1.0)

              cl_uw(J,K) = f1*(CL*dlxp(2,J)/dt)+f2*cl_uw(J,K)

              obc_uw(J,K)=uoo_wb0(0,J,K)+ cl_uw(J,K)*dtdxpo(1,J)*d01
!             if(cl_uw(J,K) .lt. obcerr)then
!                obc_uw(J,K)=uoo_wb0(1,J,K)
!             else
!                obc_uw(J,K)=uoo_wb0(0,J,K)+ cl_uw(J,K)*dtdxpo(1,J)*d01
!             endif
            enddo
          enddo
          uoo_wb2=uoo_wb1
          uoo_wb1=uoo_wb0

          voe_wb0=VVV(1:3,:,:)
          do k=1,ke
            do j=J_start,je
              d01 = ab1*(voe_wb0(2,J,K)-voe_wb0(1,J,K)) +               &
     &              ab2*(voe_wb1(2,J,K)-voe_wb1(1,J,K))
              d12 = ab1*(voe_wb1(3,J,K)-voe_wb1(2,J,K)) +               &
     &              ab2*(voe_wb2(3,J,K)-voe_wb2(2,J,K))

              IF (abs(d12) .lt. obcerr) THEN
                 CL=0.
              ELSE
                 CL=(voe_wb0(2,J,K)-voe_wb1(2,J,K)) / d12
              ENDIF
              CL=min(max(CL,0.0),1.0)

              cl_vw(J,K) = f1*(CL*dlxu(2,J)/dt)+f2*cl_vw(J,K)

              obc_vw(J,K)=voe_wb0(1,J,K)+ cl_vw(J,K)*dtdxuo(1,J)*d01
!             if(cl_vw(J,K) .lt. obcerr)then
!                obc_vw(J,K)=voe_wb0(2,J,K)
!             else
!               obc_vw(J,K)=voe_wb0(1,J,K)+ cl_vw(J,K)*dtdxuo(1,J)*d01
!             endif
            enddo
          enddo
          voe_wb2=voe_wb1
          voe_wb1=voe_wb0
       endif
      endif
      
      ! east OB
      IF ( East .AND. have_g_ie) THEN
       if (OrlanskiEast) then       
          uoo_eb0=UUU(ie2:ie,:,:)
          do k=1,ke
            do j=1,je
              d01 = ab1*(uoo_eb0(IE1,J,K)-uoo_eb0(IE,J,K)) +            &
     &              ab2*(uoo_eb1(IE1,J,K)-uoo_eb1(IE,J,K))
              d12 = ab1*(uoo_eb1(IE2,J,K)-uoo_eb1(IE1,J,K)) +           &
     &              ab2*(uoo_eb2(IE2,J,K)-uoo_eb2(IE1,J,K))

              IF (abs(d12) .lt. obcerr) THEN
                 CL=0.
              ELSE
                 CL=(uoo_eb0(IE1,J,K)-uoo_eb1(IE1,J,K)) / d12
              ENDIF
              CL=min(max(CL,0.0),1.0)

              cl_ue(J,K) = f1*(CL*dlxp(IE1,J)/dt)+f2*cl_ue(J,K)

              obc_ue(J,K)=uoo_eb0(IE,J,K)+ cl_ue(J,K)*dtdxpo(IE,J)*d01
!              if(cl_ue(J,K) .lt. obcerr)then
!                obc_ue(J,K)=uoo_eb0(IE1,J,K)
!              else
!                obc_ue(J,K)=uoo_eb0(IE,J,K)+ cl_ue(J,K)*dtdxpo(IE,J)*d01
!              endif
            enddo
          enddo
          uoo_eb2=uoo_eb1
          uoo_eb1=uoo_eb0

          voe_eb0=VVV(ie2:ie,:,:)
          do k=1,ke
            do j=J_start,je
              d01 = ab1*(voe_eb0(IE1,J,K)-voe_eb0(IE,J,K)) +            &
     &              ab2*(voe_eb1(IE1,J,K)-voe_eb1(IE,J,K))
              d12 = ab1*(voe_eb1(IE2,J,K)-voe_eb1(IE1,J,K)) +           &
     &              ab2*(voe_eb2(IE2,J,K)-voe_eb2(IE1,J,K))

              IF (abs(d12) .lt. obcerr) THEN
                 CL=0.
              ELSE
                 CL=(voe_eb0(IE1,J,K)-voe_eb1(IE1,J,K)) / d12
              ENDIF
              CL=min(max(CL,0.0),1.0)

              cl_ve(J,K) = f1*(CL*dlxu(IE2,J)/dt)+f2*cl_ve(J,K)

              obc_ve(J,K)=voe_eb0(IE,J,K)+ cl_ve(J,K)*dtdxuo(IE1,J)*d01
!              if(cl_ve(J,K) .lt. obcerr)then
!                obc_ve(J,K)=voe_eb0(IE1,J,K)
!              else
!                obc_ve(J,K)=voe_eb0(IE,J,K)+ cl_ve(J,K)*dtdxuo(IE1,J)*d01
!              endif
            enddo
          enddo
          voe_eb2=voe_eb1
          voe_eb1=voe_eb0
       endif
      endif
      
      ! north OB
      IF ( North .AND. have_g_js) THEN 
       if (OrlanskiNorth) then        
        voe_nb0=VVV(:,0:2,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(voe_nb0(I,1,K)-voe_nb0(I,0,K)) +                 &
     &            ab2*(voe_nb1(I,1,K)-voe_nb1(I,0,K))
            d12 = ab1*(voe_nb1(I,2,K)-voe_nb1(I,1,K)) +                 &
     &            ab2*(voe_nb2(I,2,K)-voe_nb2(I,1,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(voe_nb0(I,1,K)-voe_nb1(I,1,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_vn(I,K) = f1*(CL*dlyp(I,2)/dt)+f2*cl_vn(I,K)

            obc_vn(I,K)=voe_nb0(I,0,K)+ cl_vn(I,K)*dtdyo(I,1)*d01
!            if(cl_vn(I,K) .lt. obcerr)then
!              obc_vn(I,K)=voe_nb0(I,1,K)
!            else
!              obc_vn(I,K)=voe_nb0(I,0,K)+ cl_vn(I,K)*dtdyo(I,1)*d01
!            endif
          enddo
        enddo
        voe_nb2=voe_nb1
        voe_nb1=voe_nb0

        uoo_nb0=UUU(:,1:3,:)
        do k=1,ke
          do i=I_start,ie
            d01 = ab1*(uoo_nb0(I,2,K)-uoo_nb0(I,1,K)) +                 &
     &            ab2*(uoo_nb1(I,2,K)-uoo_nb1(I,1,K))
            d12 = ab1*(uoo_nb1(I,3,K)-uoo_nb1(I,2,K)) +                 &
     &            ab2*(uoo_nb2(I,3,K)-uoo_nb2(I,2,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(uoo_nb0(I,2,K)-uoo_nb1(I,2,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_un(I,K) = f1*(CL*dlyv(I,2)/dt)+f2*cl_un(I,K)

            obc_un(I,K)=uoo_nb0(I,1,K)+ cl_un(I,K)*dpye(I,1)*d01
!            if(cl_un(I,K) .lt. obcerr)then
!              obc_un(I,K)=uoo_nb0(I,2,K)
!            else
!              obc_un(I,K)=uoo_nb0(I,1,K)+ cl_un(I,K)*dpye(I,1)*d01
!            endif
          enddo
        enddo
        uoo_nb2=uoo_nb1
        uoo_nb1=uoo_nb0        
       endif
      endif
      
      ! south OB
      IF ( South .AND. have_g_je) THEN   
       if (OrlanskiSouth) then        	    
        voe_sb0=VVV(:,je2:je,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(voe_sb0(I,JE1,K)-voe_sb0(I,JE,K))  +             &
     &            ab2*(voe_sb1(I,JE1,K)-voe_sb1(I,JE,K))
            d12 = ab1*(voe_sb1(I,JE2,K)-voe_sb1(I,JE1,K)) +             &
     &            ab2*(voe_sb2(I,JE2,K)-voe_sb2(I,JE1,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(voe_sb0(I,JE1,K)-voe_sb1(I,JE1,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_vs(I,K) = f1*(CL*dlyp(I,JE1)/dt)+f2*cl_vs(I,K)

            obc_vs(I,K)=voe_sb0(I,JE,K)+ cl_vs(I,K)*dtdyo(I,JE)*d01
!            if(cl_vs(I,K) .lt. obcerr)then
!              obc_vs(I,K)=voe_sb0(I,JE1,K)
!            else
!              obc_vs(I,K)=voe_sb0(I,JE,K)+ cl_vs(I,K)*dtdyo(I,JE)*d01
!            endif
          enddo
        enddo
        voe_sb2=voe_sb1
        voe_sb1=voe_sb0

        uoo_sb0=UUU(:,je2:je,:)
        do k=1,ke
          do i=I_start,ie
            d01 = ab1*(uoo_sb0(I,JE1,K)-uoo_sb0(I,JE,K))  +             &
     &            ab2*(uoo_sb1(I,JE1,K)-uoo_sb1(I,JE,K))
            d12 = ab1*(uoo_sb1(I,JE2,K)-uoo_sb1(I,JE1,K)) +             &
     &            ab2*(uoo_sb2(I,JE2,K)-uoo_sb2(I,JE1,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(uoo_sb0(I,JE1,K)-uoo_sb1(I,JE1,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_us(I,K) = f1*(CL*dlyv(I,JE2)/dt)+f2*cl_us(I,K)

            obc_us(I,K)=uoo_sb0(I,JE,K)+ cl_us(I,K)*dpye(I,JE1)*d01
!            if(cl_us(I,K) .lt. obcerr)then
!              obc_us(I,K)=uoo_sb0(I,JE1,K)
!            else
!              obc_us(I,K)=uoo_sb0(I,JE,K)+ cl_us(I,K)*dpye(I,JE1)*d01
!            endif
          enddo
        enddo
        uoo_sb2=uoo_sb1
        uoo_sb1=uoo_sb0
       endif
      endif
      
!      CALL UPDATE_OBCS_UV(UUU,VVV)

!      UUU0(:,:,:)=UUU(:,:,:)
!      VVV0(:,:,:)=VVV(:,:,:)

      RETURN
      END SUBROUTINE PRED_OBCS_UV

 
! OBC TS forcing
      SUBROUTINE FORC_OBCS_TS(i_step)
      
      REAL CLI_TWG(JE_G,KE),CLI_TEG(JE_G,KE),CLI_TNG(IE_G,KE),CLI_TSG(IE_G,KE)
      REAL CLI_SWG(JE_G,KE),CLI_SEG(JE_G,KE),CLI_SNG(IE_G,KE),CLI_SSG(IE_G,KE)

      IF ( MOD( i_step,NINT(CLI_STEP/DT) ) .EQ.  1 ) THEN

        IF(p_pe==p_io) THEN
          DO K=1,KE
          	
          	if ( West ) then 
               READ(IO_IN_OBCTW) CLI_TWG(:,K)
               READ(IO_IN_OBCSW) CLI_SWG(:,K)
            endif
            if ( East ) then             
               READ(IO_IN_OBCTE) CLI_TEG(:,K)
               READ(IO_IN_OBCSE) CLI_SEG(:,K)
            endif  
            if ( North ) then         
               READ(IO_IN_OBCTN) CLI_TNG(:,K)
               READ(IO_IN_OBCSN) CLI_SNG(:,K)
            endif
            if ( South ) then            
               READ(IO_IN_OBCTS) CLI_TSG(:,K)
               READ(IO_IN_OBCSS) CLI_SSG(:,K)
            endif
            
          ENDDO
        ENDIF      
      

        if ( West ) then 
        	CALL p_bcast(CLI_TWG,p_io)
          CALL p_bcast(CLI_SWG,p_io)
          if (have_g_is) then
           cli_tw(:,:)=cli_twg(p_joff+1:p_joff+je,:)*weto(1,:,:)
           cli_sw(:,:)=cli_swg(p_joff+1:p_joff+je,:)*weto(1,:,:)
           obc_tw(:,:)=cli_tw(:,:)
           obc_sw(:,:)=cli_sw(:,:)
          endif
        endif
        
        if ( East ) then 
         CALL p_bcast(CLI_TEG,p_io)
         CALL p_bcast(CLI_SEG,p_io)        	         
         if (have_g_ie) then
          cli_te(:,:)=cli_teg(p_joff+1:p_joff+je,:)*weto(ie,:,:)
          cli_se(:,:)=cli_seg(p_joff+1:p_joff+je,:)*weto(ie,:,:)
          obc_te(:,:)=cli_te(:,:)
          obc_se(:,:)=cli_se(:,:)          
         endif
        endif
 
       ! if ( East .AND. have_g_ie ) then        
        ! test for NaN in occlit.f90
       !  write(IO_STDOUT,*) 'FORC_OBCS_TS'
       !  write(IO_STDOUT,*) 'obc_te(k=1)', obc_te(:,1)
       !  write(IO_STDOUT,*) 'obc_se(k=ke)', obc_se(:,ke)           
       ! endif
      
        if ( North ) then   
         CALL p_bcast(CLI_TNG,p_io)               
         CALL p_bcast(CLI_SNG,p_io)          	             
         if (have_g_js) then
          cli_tn(:,:)=cli_tng(p_ioff+1:p_ioff+ie,:)*weto(:,1,:)
          cli_sn(:,:)=cli_sng(p_ioff+1:p_ioff+ie,:)*weto(:,1,:)
          obc_tn(:,:)=cli_tn(:,:)
          obc_sn(:,:)=cli_sn(:,:)          
         endif
        endif      
        
        if ( South ) then  
         CALL p_bcast(CLI_TSG,p_io)
         CALL p_bcast(CLI_SSG,p_io)       	      
         if (have_g_je) then
          cli_ts(:,:)=cli_tsg(p_ioff+1:p_ioff+ie,:)*weto(:,je,:)
          cli_ss(:,:)=cli_ssg(p_ioff+1:p_ioff+ie,:)*weto(:,je,:)
          obc_ts(:,:)=cli_ts(:,:)
          obc_ss(:,:)=cli_ss(:,:)           
         endif
        endif
       
      ENDIF
      
      END SUBROUTINE FORC_OBCS_TS

! relax OBC T&S forcing to model fields (with sponge layer)
      SUBROUTINE SPONGE_OBCS_TS(TTT0,SSS0)
                 
      REAL TTT(IE,JE,KE),TTT0(IE,JE,KE),SSS(IE,JE,KE),SSS0(IE,JE,KE)
      REAL trelaxw(JE),trelaxe(JE),trelaxn(IE),trelaxs(IE)
      REAL srelaxw(JE),srelaxe(JE),srelaxn(IE),srelaxs(IE)      
      REAL lambda     

      TTT(:,:,:)=TTT0(:,:,:)
      SSS(:,:,:)=SSS0(:,:,:)      


      TFKO=0.0
      TFKO_W=0.0
      TFKO_E=0.0
      TFKO_N=0.0
      TFKO_S=0.0
      
      SFKO=0.0
      SFKO_W=0.0
      SFKO_E=0.0
      SFKO_N=0.0
      SFKO_S=0.0      

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
        DO k=1,ke      	
        DO i=1,spongew
          s1=real(spongew-i+1)
          s2=real(i-1)
          trelaxw(:)=(s1*obc_tw(:,k)+s2*TTT(i,:,k))/real(spongew)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongew)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          TFKO_W(i,:,k)=weto(1,:,k)*weto(i,:,k)*lambda*(trelaxw(:)-TTT(i,:,k))
        ENDDO
        
        DO i=1,spongew
          s1=real(spongew-i+1)
          s2=real(i-1)
          srelaxw(:)=(s1*obc_sw(:,k)+s2*SSS(i,:,k))/real(spongew)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongew)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          SFKO_W(i,:,k)=weto(1,:,k)*weto(i,:,k)*lambda*(srelaxw(:)-SSS(i,:,k))
        ENDDO        
        ENDDO       
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
       	
       	DO k=1,ke
        DO i=ie,ie-spongee+1,-1
          s1=real(i-ie+spongee)
          s2=real(ie-i)
          trelaxe(:)=(s1*obc_te(:,k)+s2*TTT(i,:,k))/real(spongee)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongee)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          TFKO_E(i,:,k)=weto(ie,:,k)*weto(i,:,k)*lambda*(trelaxe(:)-TTT(i,:,k))
        ENDDO
 
        DO i=ie,ie-spongee+1,-1
          s1=real(i-ie+spongee)
          s2=real(ie-i)
          srelaxe(:)=(s1*obc_se(:,k)+s2*SSS(i,:,k))/real(spongee)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongee)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          SFKO_E(i,:,k)=weto(ie,:,k)*weto(i,:,k)*lambda*(srelaxe(:)-SSS(i,:,k))
        ENDDO      
      
        ENDDO     
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
       	DO k=1,ke
        DO j=1,spongen
          s1=real(spongen-j+1)
          s2=real(j-1)
          trelaxn(:)=(s1*obc_tn(:,k)+s2*TTT(:,j,k))/real(spongen)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongen)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          TFKO_N(:,j,k)=weto(:,1,k)*weto(:,j,k)*lambda*(trelaxn(:)-TTT(:,j,k))
        ENDDO

        DO j=1,spongen
          s1=real(spongen-j+1)
          s2=real(j-1)
          srelaxn(:)=(s1*obc_sn(:,k)+s2*SSS(:,j,k))/real(spongen)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(spongen)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          SFKO_N(:,j,k)=weto(:,1,k)*weto(:,j,k)*lambda*(srelaxn(:)-SSS(:,j,k))
        ENDDO
        
        ENDDO        
       endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
       	DO k=1,ke       	
        DO j=je,je-sponges+1,-1
          s1=real(j-je+sponges)
          s2=real(je-j)
          trelaxs(:)=(s1*obc_ts(:,k)+s2*TTT(:,j,k))/real(sponges)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(sponges)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          TFKO_S(:,j,k)=weto(:,je,k)*weto(:,j,k)*lambda*(trelaxs(:)-TTT(:,j,k))
        ENDDO

        DO j=je,je-sponges+1,-1
          s1=real(j-je+sponges)
          s2=real(je-j)
          srelaxs(:)=(s1*obc_ss(:,k)+s2*SSS(:,j,k))/real(sponges)
          lambda=(s1*rel_tsb+s2*rel_tsi)/real(sponges)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          SFKO_S(:,j,k)=weto(:,je,k)*weto(:,j,k)*lambda*(srelaxs(:)-SSS(:,j,k))
        ENDDO
        ENDDO        
        
       endif
      endif

      TFKO=TFKO_W+TFKO_E+TFKO_N+TFKO_S
      SFKO=SFKO_W+SFKO_E+SFKO_N+SFKO_S      
      

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then       	 
        DO i=1,spongew
          DO j=1,spongen
            if (i<j) then
              TFKO(i,j,:)=TFKO_W(i,j,:)
              SFKO(i,j,:)=SFKO_W(i,j,:)             
            else
              TFKO(i,j,:)=TFKO_N(i,j,:)
              SFKO(i,j,:)=SFKO_N(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then        	
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              TFKO(i,j,:)=TFKO_W(i,j,:)
              SFKO(i,j,:)=SFKO_W(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_S(i,j,:)
              SFKO(i,j,:)=SFKO_S(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then     
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              TFKO(i,j,:)=TFKO_E(i,j,:)
              SFKO(i,j,:)=SFKO_E(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_N(i,j,:)
              SFKO(i,j,:)=SFKO_N(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then        	
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              TFKO(i,j,:)=TFKO_E(i,j,:)
              SFKO(i,j,:)=SFKO_E(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_S(i,j,:)
              SFKO(i,j,:)=SFKO_S(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      TTT0=TTT+TFKO
      SSS0=SSS+SFKO      

      CALL bounds_exch(TTT0)
      CALL bounds_exch(SSS0)
      
      ! if ( East .AND. have_g_ie ) then          
      !  write(IO_STDOUT,*) 'SPONGE_OBCS_TS'
      !  write(IO_STDOUT,*) 'TTT0(k=1)', TTT0(IE,:,1)
      !  write(IO_STDOUT,*) 'TTT0(k=ke)', TTT0(IE,:,ke)  
      !  write(IO_STDOUT,*) 'SSS0(k=1)', SSS0(IE,:,1)
      !  write(IO_STDOUT,*) 'SSS0(k=ke)', SSS0(IE,:,ke)                      
      ! ENDIF

      END SUBROUTINE SPONGE_OBCS_TS

      
! Open Boundary Radiation Condition ts
      SUBROUTINE PRED_OBCS_TS(TTT0,SSS0)

      REAL TTT0(IE,JE,KE),TTT(IE,JE,KE)
      REAL SSS0(IE,JE,KE),SSS(IE,JE,KE)

      TTT(:,:,:)=TTT0(:,:,:)
      SSS(:,:,:)=SSS0(:,:,:)
      
!  Added code to allow filtering of phase speed following
!  Blumberg and Kantha. There is now a separate array
!  cl_**, where **=Variable(U,V,T,S,W)Boundary(E,W,N,S) for
!  the dimensional phase speed. These arrays are initialized to
!  zero. CVEL_** is filtered according to
!  cl_** = f1*cl_**(new)+(1-f1)*cl_**(old).
!  f1=1.0 turns off filtering.
      f1 = dt/time_scale
      f2 = 1.0-f1
!  Changed code to average phase speed. A new variable
!  'time_scale' was created. This variable must now be
!  specified. Then, f1=dT/time_scale.Since the goal is 
!  to smooth out the 'singularities' in thediagnosed phase 
!  speed, time_scale could be picked as the duration of the
!  singular period in the unfiltered case. Thus,for a plane wave
!  time_scale might be the time take for thewave to travel 
!  a distance DX, where DX is the width of the region near which 
!  d(phi)/dx is small.  

      ! west OB
      if ( West ) then      
       if (have_g_is) then
        tho_wb0=TTT(1:3,:,:)
        do k=1,ke
          do j=1,je
            d01 = ab1*(tho_wb0(2,J,K)-tho_wb0(1,J,K)) +                 &
     &            ab2*(tho_wb1(2,J,K)-tho_wb1(1,J,K))
            d12 = ab1*(tho_wb1(3,J,K)-tho_wb1(2,J,K)) +                 &
     &            ab2*(tho_wb2(3,J,K)-tho_wb2(2,J,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(tho_wb0(2,J,K)-tho_wb1(2,J,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_tw(J,K) = f1*(CL*dlxu(2,J)/dt)+f2*cl_tw(J,K)

!            obc_tw(J,K)=tho_wb0(1,J,K)+ cl_tw(J,K)*dtdxuo(1,J)*d01
            if(cl_tw(J,K) .lt. obcerr)then
              obc_tw(J,K)=tho_wb0(1,J,K)+                               &
     &                       dt*obrel_tem*(cli_tw(J,K)-tho_wb0(1,J,K))
            else
              obc_tw(J,K)=tho_wb0(1,J,K)+ cl_tw(J,K)*dtdxuo(1,J)*d01
            endif
          enddo
        enddo
        tho_wb2=tho_wb1
        tho_wb1=tho_wb0

        sao_wb0=SSS(1:3,:,:)
        do k=1,ke
          do j=1,je
            d01 = ab1*(sao_wb0(2,J,K)-sao_wb0(1,J,K)) +                 &
     &            ab2*(sao_wb1(2,J,K)-sao_wb1(1,J,K))
            d12 = ab1*(sao_wb1(3,J,K)-sao_wb1(2,J,K)) +                 &
     &            ab2*(sao_wb2(3,J,K)-sao_wb2(2,J,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(sao_wb0(2,J,K)-sao_wb1(2,J,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_sw(J,K) = f1*(CL*dlxu(2,J)/dt)+f2*cl_sw(J,K)

!            obc_sw(J,K)=sao_wb0(1,J,K)+ cl_sw(J,K)*dtdxuo(1,J)*d01
            if(cl_sw(J,K) .lt. obcerr)then
              obc_sw(J,K)=sao_wb0(1,J,K)+                               &
     &                       dt*obrel_sal*(cli_sw(J,K)-sao_wb0(1,J,K))
            else
              obc_sw(J,K)=sao_wb0(1,J,K)+ cl_sw(J,K)*dtdxuo(1,J)*d01
            endif
          enddo
        enddo
        sao_wb2=sao_wb1
        sao_wb1=sao_wb0
       endif
      endif
      
      ! east OB
      if ( East ) then       
       if (have_g_ie) then
        tho_eb0=TTT(ie2:ie,:,:)
        do k=1,ke
          do j=1,je
            d01 = ab1*(tho_eb0(IE1,J,K)-tho_eb0(IE,J,K)) +              &
     &            ab2*(tho_eb1(IE1,J,K)-tho_eb1(IE,J,K))
            d12 = ab1*(tho_eb1(IE2,J,K)-tho_eb1(IE1,J,K)) +             &
     &            ab2*(tho_eb2(IE2,J,K)-tho_eb2(IE1,J,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(tho_eb0(IE1,J,K)-tho_eb1(IE1,J,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_te(J,K) = f1*(CL*dlxu(IE2,J)/dt)+f2*cl_te(J,K)

!            obc_te(J,K)=tho_eb0(IE,J,K)+ cl_te(J,K)*dtdxuo(IE1,J)*d01
            if(cl_te(J,K) .lt. obcerr)then
              obc_te(J,K)=tho_eb0(IE,J,K)+                              &
     &                       dt*obrel_tem*(cli_te(J,K)-tho_eb0(IE,J,K))
            else
              obc_te(J,K)=tho_eb0(IE,J,K)+ cl_te(J,K)*dtdxuo(IE1,J)*d01
            endif
          enddo
        enddo
        tho_eb2=tho_eb1
        tho_eb1=tho_eb0

        sao_eb0=SSS(ie2:ie,:,:)
        do k=1,ke
          do j=1,je
            d01 = ab1*(sao_eb0(IE1,J,K)-sao_eb0(IE,J,K)) +              &
     &            ab2*(sao_eb1(IE1,J,K)-sao_eb1(IE,J,K))
            d12 = ab1*(sao_eb1(IE2,J,K)-sao_eb1(IE1,J,K)) +             &
     &            ab2*(sao_eb2(IE2,J,K)-sao_eb2(IE1,J,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(sao_eb0(IE1,J,K)-sao_eb1(IE1,J,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_se(J,K) = f1*(CL*dlxu(IE2,J)/dt)+f2*cl_se(J,K)

!            obc_se(J,K)=sao_eb0(IE,J,K)+ cl_se(J,K)*dtdxuo(IE1,J)*d01
            if(cl_se(J,K) .lt. obcerr)then
              obc_se(J,K)=sao_eb0(IE,J,K)+                              &
     &                       dt*obrel_sal*(cli_se(J,K)-sao_eb0(IE,J,K))
            else
              obc_se(J,K)=sao_eb0(IE,J,K)+ cl_se(J,K)*dtdxuo(IE1,J)*d01
            endif
          enddo
        enddo
        sao_eb2=sao_eb1
        sao_eb1=sao_eb0
       endif
      endif
      
      ! north OB
      if ( North ) then      
       if (have_g_js) then
        tho_nb0=TTT(:,1:3,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(tho_nb0(I,2,K)-tho_nb0(I,1,K)) +                 &
     &            ab2*(tho_nb1(I,2,K)-tho_nb1(I,1,K))
            d12 = ab1*(tho_nb1(I,3,K)-tho_nb1(I,2,K)) +                 &
     &            ab2*(tho_nb2(I,3,K)-tho_nb2(I,2,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(tho_nb0(I,2,K)-tho_nb1(I,2,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_tn(I,K) = f1*(CL*dlyv(I,2)/dt)+f2*cl_tn(I,K)

!            obc_tn(I,K)=tho_nb0(I,1,K)+ cl_tn(I,K)*dpye(I,1)*d01
            if(cl_tn(I,K) .lt. obcerr)then
              obc_tn(I,K)=tho_nb0(I,1,K)+                               &
     &                       dt*obrel_tem*(cli_tn(I,K)-tho_nb0(I,1,K))
            else
              obc_tn(I,K)=tho_nb0(I,1,K)+ cl_tn(I,K)*dpye(I,1)*d01
            endif
          enddo
        enddo
        tho_nb2=tho_nb1
        tho_nb1=tho_nb0

        sao_nb0=SSS(:,1:3,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(sao_nb0(I,2,K)-sao_nb0(I,1,K)) +                 &
     &            ab2*(sao_nb1(I,2,K)-sao_nb1(I,1,K))
            d12 = ab1*(sao_nb1(I,3,K)-sao_nb1(I,2,K)) +                 &
     &            ab2*(sao_nb2(I,3,K)-sao_nb2(I,2,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(sao_nb0(I,2,K)-sao_nb1(I,2,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_sn(I,K) = f1*(CL*dlyv(I,2)/dt)+f2*cl_sn(I,K)

!            obc_sn(I,K)=sao_nb0(I,1,K)+ cl_sn(I,K)*dpye(I,1)*d01
            if(cl_sn(I,K) .lt. obcerr)then
              obc_sn(I,K)=sao_nb0(I,1,K)+                               &
     &                       dt*obrel_sal*(cli_sn(I,K)-sao_nb0(I,1,K))
            else
              obc_sn(I,K)=sao_nb0(I,1,K)+ cl_sn(I,K)*dpye(I,1)*d01
            endif
          enddo
        enddo
        sao_nb2=sao_nb1
        sao_nb1=sao_nb0
       endif
      endif
      
      ! south OB
      if ( South ) then       
       if (have_g_je) then
        tho_sb0=TTT(:,je2:je,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(tho_sb0(I,JE1,K)-tho_sb0(I,JE,K))  +             &
     &            ab2*(tho_sb1(I,JE1,K)-tho_sb1(I,JE,K))
            d12 = ab1*(tho_sb1(I,JE2,K)-tho_sb1(I,JE1,K)) +             &
     &            ab2*(tho_sb2(I,JE2,K)-tho_sb2(I,JE1,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(tho_sb0(I,JE1,K)-tho_sb1(I,JE1,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_ts(I,K) = f1*(CL*dlyv(I,JE2)/dt)+f2*cl_ts(I,K)

!            obc_ts(I,K)=tho_sb0(I,JE,K)+ cl_ts(I,K)*dpye(I,JE1)*d01
            if(cl_ts(I,K) .lt. obcerr)then
              obc_ts(I,K)=tho_sb0(I,JE,K)+                              &
     &                       dt*obrel_tem*(cli_ts(I,K)-tho_sb0(I,JE,K))
            else
              obc_ts(I,K)=tho_sb0(I,JE,K)+ cl_ts(I,K)*dpye(I,JE1)*d01
            endif
          enddo
        enddo
        tho_sb2=tho_sb1
        tho_sb1=tho_sb0

        sao_sb0=SSS(:,je2:je,:)
        do k=1,ke
          do i=1,ie
            d01 = ab1*(sao_sb0(I,JE1,K)-sao_sb0(I,JE,K))  +             &
     &            ab2*(sao_sb1(I,JE1,K)-sao_sb1(I,JE,K))
            d12 = ab1*(sao_sb1(I,JE2,K)-sao_sb1(I,JE1,K)) +             &
     &            ab2*(sao_sb2(I,JE2,K)-sao_sb2(I,JE1,K))

            IF (abs(d12) .lt. obcerr) THEN
               CL=0.
            ELSE
               CL=(sao_sb0(I,JE1,K)-sao_sb1(I,JE1,K)) / d12
            ENDIF
            CL=min(max(CL,0.0),1.0)

            cl_ss(I,K) = f1*(CL*dlyv(I,JE2)/dt)+f2*cl_ss(I,K)

!            obc_ss(I,K)=sao_sb0(I,JE,K)+ cl_ss(I,K)*dpye(I,JE1)*d01
            if(cl_ss(I,K) .lt. obcerr)then
              obc_ss(I,K)=sao_sb0(I,JE,K)+                              &
     &                       dt*obrel_sal*(cli_ss(I,K)-sao_sb0(I,JE,K))
            else
              obc_ss(I,K)=sao_sb0(I,JE,K)+ cl_ss(I,K)*dpye(I,JE1)*d01
            endif
          enddo
        enddo
        sao_sb2=sao_sb1
        sao_sb1=sao_sb0
       endif
      endif
      
      CALL UPDATE_OBCS_TS

      END SUBROUTINE PRED_OBCS_TS


! relax OBC eta forcing to model fields (with sponge layer)
      SUBROUTINE SPONGE_OBCS_ETA(ZZZ0)

      REAL zrelaxw(JE),zrelaxe(JE),zrelaxn(IE),zrelaxs(IE)
      REAL lambda,ZZZ0(IE,JE),ZZZ(IE,JE)
      
      ZZZ(:,:) = ZZZ0(:,:)

      ZFO=0.0
      ZFO_W=0.0
      ZFO_E=0.0
      ZFO_N=0.0
      ZFO_S=0.0

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
        DO i=1,spongew
          s1=real(spongew-i+1)
          s2=real(i-1)
          zrelaxw(:)=(s1*obc_zw(:)+s2*ZZZ(i,:))/real(spongew)
          lambda=(s1*rel_etab+s2*rel_etai)/real(spongew)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          ZFO_W(i,:)=weto(1,:,1)*weto(i,:,1)*lambda*(zrelaxw(:)-ZZZ(i,:))
        ENDDO
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
        DO i=ie,ie-spongee+1,-1
          s1=real(i-ie+spongee)
          s2=real(ie-i)
          zrelaxe(:)=(s1*obc_ze(:)+s2*ZZZ(i,:))/real(spongee)
          lambda=(s1*rel_etab+s2*rel_etai)/real(spongee)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          ZFO_E(i,:)=weto(ie,:,1)*weto(i,:,1)*lambda*(zrelaxe(:)-ZZZ(i,:))
        ENDDO
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
        DO j=1,spongen
          s1=real(spongen-j+1)
          s2=real(j-1)
          zrelaxn(:)=(s1*obc_zn(:)+s2*ZZZ(:,j))/real(spongen)
          lambda=(s1*rel_etab+s2*rel_etai)/real(spongen)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          ZFO_N(:,j)=weto(:,1,1)*weto(:,j,1)*lambda*(zrelaxn(:)-ZZZ(:,j))
        ENDDO
       endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
        DO j=je,je-sponges+1,-1
          s1=real(j-je+sponges)
          s2=real(je-j)
          zrelaxs(:)=(s1*obc_zs(:)+s2*ZZZ(:,j))/real(sponges)
          lambda=(s1*rel_etab+s2*rel_etai)/real(sponges)
          IF (lambda.ne.0.) THEN
            lambda=dt/lambda
          ELSE
            lambda=0.0
          ENDIF
          ZFO_S(:,j)=weto(:,je,1)*weto(:,j,1)*lambda*(zrelaxs(:)-ZZZ(:,j))
        ENDDO
       endif
      endif

      ZFO=ZFO_W+ZFO_E+ZFO_N+ZFO_S

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then       	 
        DO i=1,spongew
          DO j=1,spongen
            if (i<j) then
              ZFO(i,j)=ZFO_W(i,j)
            else
              ZFO(i,j)=ZFO_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then        	
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              ZFO(i,j)=ZFO_W(i,j)
            else
              ZFO(i,j)=ZFO_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then     
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              ZFO(i,j)=ZFO_E(i,j)
            else
              ZFO(i,j)=ZFO_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then        	
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              ZFO(i,j)=ZFO_E(i,j)
            else
              ZFO(i,j)=ZFO_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ZZZ0=ZZZ+ZFO

      CALL bounds_exch(ZZZ0)

      END SUBROUTINE SPONGE_OBCS_ETA


! update OBC predict uv to model fields
      SUBROUTINE UPDATE_OBCS_UV(UUU0,VVV0)

      REAL UUU0(1:IE-I_start+1,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,1:JE-J_start+1,KE),VVV(IE,J_start:JE,KE)

      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)

      ! west OB
      if ( West .AND. have_g_is) then        
        UUU(0,:,:)=obc_uw(:,:)*amsuo(0,:,:)
        VVV(1,:,:)=obc_vw(:,:)*amsue(1,:,:)
      endif

      ! east OB
      if ( East .AND. have_g_ie) then       
        UUU(ie,:,:)=obc_ue(:,:)*amsuo(ie,:,:)
        VVV(ie,:,:)=obc_ve(:,:)*amsue(ie,:,:)
      endif

      ! north OB
      if ( North .AND. have_g_js) then
        UUU(:,1,:)=obc_un(:,:)*amsuo(:,1,:)
        VVV(:,0,:)=obc_vn(:,:)*amsue(:,0,:)
      endif

      ! south OB
      if ( South .AND. have_g_je) then     
        UUU(:,je,:)=obc_us(:,:)*amsuo(:,je,:)
        VVV(:,je,:)=obc_vs(:,:)*amsue(:,je,:)
      endif

      CALL bounds_exch(UUU)
      CALL bounds_exch(VVV)

      UUU0(:,:,:)=UUU(:,:,:)
      VVV0(:,:,:)=VVV(:,:,:)

      RETURN
      END SUBROUTINE UPDATE_OBCS_UV

      SUBROUTINE UPDATE_OBCS_UV_2D(UUU,VVV)
      
      REAL UUU0(1:IE-I_start+1,JE),UUU(I_start:IE,JE)
      REAL VVV0(IE,1:JE-J_start+1),VVV(IE,J_start:JE)

      UUU(:,:)=UUU0(:,:)
      VVV(:,:)=VVV0(:,:)

      ! west OB
      if (WEST) then
      if (have_g_is .AND. .NOT. SpongeWest) then        
        UUU(0,:)=obc_uzw(:)*amsuo(0,:,1)
        VVV(1,:)=obc_vzw(:)*amsue(1,:,1)
      endif
      endif

      ! east OB
      if (EAST) then
      if (have_g_ie .AND. .NOT. SpongeEast) then       
        UUU(ie,:)=obc_uze(:)*amsuo(ie,:,1)
        VVV(ie,:)=obc_vze(:)*amsue(ie,:,1)
      endif
      endif

      ! north OB
      if (NORTH) then
      if (have_g_js .AND. .NOT. SpongeNorth) then
        UUU(:,1)=obc_uzn(:)*amsuo(:,1,1)
        VVV(:,0)=obc_vzn(:)*amsue(:,0,1)
      endif
      ENDIF

      ! south OB
      if (SOUTH) then
      if ( have_g_je .AND. .NOT. SpongeSouth) then     
        UUU(:,je)=obc_uzs(:)*amsuo(:,je,1)
        VVV(:,je)=obc_vzs(:)*amsue(:,je,1)
      endif
      ENDIF

      CALL bounds_exch(UUU)
      CALL bounds_exch(VVV)

      UUU0(:,:)=UUU(:,:)
      VVV0(:,:)=VVV(:,:)
         
      RETURN
      END SUBROUTINE UPDATE_OBCS_UV_2D

! update OBC predict/relax ts to model fields
      SUBROUTINE UPDATE_OBCS_TS

      ! west OB
      if ( West ) then 
      if (have_g_is) then
        do k=1,ke
          do j=1,je
            if (weto(1,j,k).eq.1)then
              tho(1,j,k)=obc_tw(j,k)
              sao(1,j,k)=obc_sw(j,k)
            endif
          enddo
        enddo
      endif
      endif

      ! east OB
      if ( East ) then       
      if (have_g_ie) then
        do k=1,ke
          do j=1,je
            if (weto(ie,j,k).eq.1)then
              tho(ie,j,k)=obc_te(j,k)
              sao(ie,j,k)=obc_se(j,k)
            endif
          enddo
        enddo
      endif
      endif

      ! north OB
      if ( North ) then       
      if (have_g_js) then
        do k=1,ke
          do i=1,ie
            if (weto(i,1,k).eq.1)then
              tho(i,1,k)=obc_tn(i,k)
              sao(i,1,k)=obc_sn(i,k)
            endif
          enddo
        enddo
      endif
      endif

      ! south OB
      if ( South ) then 
      if (have_g_je) then
        do k=1,ke
          do i=1,ie
            if (weto(i,je,k).eq.1)then
              tho(i,je,k)=obc_ts(i,k)
              sao(i,je,k)=obc_ss(i,k)
            endif
          enddo
        enddo
      endif
      endif

      CALL bounds_exch(THO)
      CALL bounds_exch(SAO)

      END SUBROUTINE UPDATE_OBCS_TS


! CALL FORC_OBCS_UV((LDTRUN*DT))
! OBC uv forcing---Forced by velocity UV (m/s)
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
        if (add_tide) then        
         obc_uzw(:) = obc_uzw(:)*DEUTO(0,1:je)
         obc_vzw(:) = obc_vzw(:)*DEUTE(1,J_start:je)
       endif
      endif
     
     ! if ( west .AND. have_g_is ) then        
     !    write(IO_STDOUT,*) 'FORC_OBCS_UV'
     !    write(IO_STDOUT,*) 'LDTRUN', i_step
     !    write(IO_STDOUT,*) 'obc_uzw', obc_uzw(:)
     !    write(IO_STDOUT,*) 'obc_vzw', obc_vzw(:)       
     ! endif

 
      
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
        if (add_tide) then
         obc_uze(:) = obc_uze(:)*DEUTO(ie,1:je)
         obc_vze(:) = obc_vze(:)*DEUTE(ie,J_start:je)   
        endif    
      endif

     ! if ( East .AND. have_g_ie ) then
     !    write(IO_STDOUT,*) 'FORC_OBCS_UV'
     !    write(IO_STDOUT,*) 'LDTRUN', i_step
     !    write(IO_STDOUT,*) 'obc_uze', obc_uze(:)
     !    write(IO_STDOUT,*) 'obc_vze', obc_vze(:)
     ! endif


      
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
        if (add_tide) then
         obc_uzn(:) = obc_uzn(:)*DEUTO(I_start:ie,1)
         obc_vzn(:) = obc_vzn(:)*DEUTE(1:ie,0)   
        endif       
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
        if (add_tide) then
          obc_uzs(:) = obc_uzs(:)*DEUTO(I_start:ie,je)
          obc_vzs(:) = obc_vzs(:)*DEUTE(1:ie,je) 
        endif
      endif


! H.H --For total velocity Nudging
       CALL UPGRADE_OBCS_UV   
#endif
       
! H.H --For barotropic velocity Nudging
      ! if (add_tide) then
      ! CALL SPONGE_OBCS_UV_2D
      ! endif

      END SUBROUTINE FORC_OBCS_UV


      SUBROUTINE UPGRADE_OBCS_UV
          
      ! west OB
      if ( West .AND. have_g_is) then   
       do k = 1,ke
         obc_uw(:,k)=obc_uzw(:)*amsuo(0,:,k)
         obc_vw(:,k)=obc_vzw(:)*amsue(1,:,k)
        enddo
      endif
      
      ! east OB
      if ( East .AND. have_g_ie ) then  
       do k = 1,ke
         obc_ue(:,k)=obc_uze(:)*amsuo(ie,:,k)
         obc_ve(:,k)=obc_vze(:)*amsue(ie,:,k)   
        ENDDO
      ENDIF
      

      ! north OB
      if ( North .AND. have_g_js) then
       do k = 1,ke
        obc_un(:,k)=obc_uzn(:)*amsuo(:,1,k)
        obc_vn(:,k)=obc_vzn(:)*amsue(:,0,k) 
       enddo
      endif

      ! south OB
      if ( South .AND. have_g_je) then 
       do k = 1,ke      	
        obc_us(:,k)=obc_uzs(:)*amsuo(:,je,k)
        obc_vs(:,k)=obc_vzs(:)*amsue(:,je,k)
       enddo        
      endif            

      END SUBROUTINE UPGRADE_OBCS_UV
      
      
! relax OBC uv(2D) forcing to model fields (with sponge layer)
      SUBROUTINE SPONGE_OBCS_UV_2D

      REAL urelaxw(JE),urelaxe(JE),urelaxn(I_start:IE),urelaxs(I_start:IE)
      REAL vrelaxw(J_start:JE),vrelaxe(J_start:JE),vrelaxn(IE),vrelaxs(IE)
      REAL lambda

      UFO=0.0
      UFO_W=0.0
      UFO_E=0.0
      UFO_N=0.0
      UFO_S=0.0

      VFE=0.0
      VFE_W=0.0
      VFE_E=0.0
      VFE_N=0.0
      VFE_S=0.0

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
          DO i=0,spongew-1
            s1=real(spongew-i)
            s2=real(i)
            urelaxw(:)=(s1*obc_uzw(:)+s2*UZO(i,:))/real(spongew)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongew)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFO_W(i,:)=amsuo(0,:,1)*amsuo(i,:,1)*lambda*(urelaxw(:)-UZO(i,:))
          ENDDO

          DO i=1,spongew
            s1=real(spongew-i+1)
            s2=real(i-1)
            vrelaxw(:)=(s1*obc_vzw(:)+s2*VZE(i,:))/real(spongew)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongew)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFE_W(i,:)=amsue(1,:,1)*amsue(i,:,1)*lambda*(vrelaxw(:)-VZE(i,:))
          ENDDO
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
          DO i=ie,ie-spongee+1,-1
            s1=real(i-ie+spongee)
            s2=real(ie-i)
            urelaxe(:)=(s1*obc_uze(:)+s2*UZO(i,:))/real(spongee)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongee)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFO_E(i,:)=amsuo(ie,:,1)*amsuo(i,:,1)*lambda*(urelaxe(:)-UZO(i,:))
          ENDDO

          DO i=ie,ie-spongee+1,-1
            s1=real(i-ie+spongee)
            s2=real(ie-i)
            vrelaxe(:)=(s1*obc_vze(:)+s2*VZE(i,:))/real(spongee)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongee)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFE_E(i,:)=amsue(ie,:,1)*amsue(i,:,1)*lambda*(vrelaxe(:)-VZE(i,:))
          ENDDO
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
          DO j=1,spongen
            s1=real(spongen-j+1)
            s2=real(j-1)
            urelaxn(:)=(s1*obc_uzn(:)+s2*UZO(:,j))/real(spongen)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongen)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFO_N(:,j)=amsuo(:,1,1)*amsuo(:,j,1)*lambda*(urelaxn(:)-UZO(:,j))
          ENDDO

          DO j=0,spongen-1
            s1=real(spongen-j)
            s2=real(j)
            vrelaxn(:)=(s1*obc_vzn(:)+s2*VZE(:,j))/real(spongen)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongen)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFE_N(:,j)=amsue(:,0,1)*amsue(:,j,1)*lambda*(vrelaxn(:)-VZE(:,j))
          ENDDO
       endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
          DO j=je,je-sponges+1,-1
            s1=real(j-je+sponges)
            s2=real(je-j)
            urelaxs(:)=(s1*obc_uzs(:)+s2*UZO(:,j))/real(sponges)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(sponges)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFO_S(:,j)=amsuo(:,je,1)*amsuo(:,j,1)*lambda*(urelaxs(:)-UZO(:,j))
          ENDDO

          DO j=je,je-sponges+1,-1
            s1=real(j-je+sponges)
            s2=real(je-j)
            vrelaxs(:)=(s1*obc_vzs(:)+s2*VZE(:,j))/real(sponges)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(sponges)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFE_S(:,j)=amsue(:,je,1)*amsue(:,j,1)*lambda*(vrelaxs(:)-VZE(:,j))
          ENDDO
       endif
      endif

      UFO=UFO_W+UFO_E+UFO_N+UFO_S
      VFE=VFE_W+VFE_E+VFE_N+VFE_S

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then  
        DO i=0,spongew-1
          DO j=1,spongen
            if ((i+1)<j) then
              UFO(i,j)=UFO_W(i,j)
            else
              UFO(i,j)=UFO_N(i,j)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=0,spongen-1
            if (i<(j+1)) then
              VFE(i,j)=VFE_W(i,j)
            else
              VFE(i,j)=VFE_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then 
        DO i=0,spongew-1
          DO j=je,je-sponges+1,-1
            if ((i+1)<(je-j+1)) then
              UFO(i,j)=UFO_W(i,j)
            else
              UFO(i,j)=UFO_S(i,j)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              VFE(i,j)=VFE_W(i,j)
            else
              VFE(i,j)=VFE_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then 
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              UFO(i,j)=UFO_E(i,j)
            else
              UFO(i,j)=UFO_N(i,j)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=0,spongen-1
            if ((ie-i+1)<(j+1)) then
              VFE(i,j)=VFE_E(i,j)
            else
              VFE(i,j)=VFE_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then   
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              UFO(i,j)=UFO_E(i,j)
            else
              UFO(i,j)=UFO_S(i,j)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              VFE(i,j)=VFE_E(i,j)
            else
              VFE(i,j)=VFE_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      UZO=UZO+UFO
      VZE=VZE+VFE

      CALL bounds_exch(UZO)
      CALL bounds_exch(VZE)

      END SUBROUTINE SPONGE_OBCS_UV_2D


! relax OBC uv(3D) forcing to model fields (with sponge layer)
      SUBROUTINE SPONGE_OBCS_UV(UUU0,VVV0)

      REAL UUU0(I_start:IE,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,J_start:JE,KE),VVV(IE,J_start:JE,KE)
      REAL urelaxw(JE),urelaxe(JE),urelaxn(I_start:IE),urelaxs(I_start:IE)
      REAL vrelaxw(J_start:JE),vrelaxe(J_start:JE),vrelaxn(IE),vrelaxs(IE)
      REAL lambda

      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)
      
      UFKO=0.0
      UFKO_W=0.0
      UFKO_E=0.0
      UFKO_N=0.0
      UFKO_S=0.0

      VFKE=0.0
      VFKE_W=0.0
      VFKE_E=0.0
      VFKE_N=0.0
      VFKE_S=0.0

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
      	  DO k=1,ke
           DO i=0,spongew-1
             s1=real(spongew-i)
             s2=real(i)
             urelaxw(:)=(s1*obc_uw(:,k)+s2*UUU(i,:,k))/real(spongew)
             lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongew)
             IF (lambda.ne.0.) THEN
               lambda=dt/lambda
             ELSE
               lambda=0.0
             ENDIF
             UFKO_W(i,:,k)=amsuo(0,:,k)*amsuo(i,:,k)*lambda*(urelaxw(:)-UUU(i,:,k))
           ENDDO
         

           DO i=1,spongew
             s1=real(spongew-i+1)
             s2=real(i-1)
             vrelaxw(:)=(s1*obc_vw(:,k)+s2*VVV(i,:,k))/real(spongew)
             lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongew)
             IF (lambda.ne.0.) THEN
               lambda=dt/lambda
             ELSE
               lambda=0.0
             ENDIF
             VFKE_W(i,:,k)=amsue(1,:,k)*amsue(i,:,k)*lambda*(vrelaxw(:)-VVV(i,:,k))
           ENDDO
          ENDDO
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
       	 DO k = 1,ke
          DO i=ie,ie-spongee+1,-1
            s1=real(i-ie+spongee)
            s2=real(ie-i)
            urelaxe(:)=(s1*obc_ue(:,k)+s2*UUU(i,:,k))/real(spongee)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongee)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFKO_E(i,:,k)=amsuo(ie,:,k)*amsuo(i,:,k)*lambda*(urelaxe(:)-UUU(i,:,k))
          ENDDO

          DO i=ie,ie-spongee+1,-1
            s1=real(i-ie+spongee)
            s2=real(ie-i)
            vrelaxe(:)=(s1*obc_ve(:,k)+s2*VVV(i,:,k))/real(spongee)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongee)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFKE_E(i,:,k)=amsue(ie,:,k)*amsue(i,:,k)*lambda*(vrelaxe(:)-VVV(i,:,k))
          ENDDO
         endDO
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
       	 DO k = 1,ke      	
          DO j=1,spongen
            s1=real(spongen-j+1)
            s2=real(j-1)
            urelaxn(:)=(s1*obc_un(:,k)+s2*UUU(:,j,k))/real(spongen)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongen)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFKO_N(:,j,k)=amsuo(:,1,k)*amsuo(:,j,k)*lambda*(urelaxn(:)-UUU(:,j,k))
          ENDDO

          DO j=0,spongen-1
            s1=real(spongen-j)
            s2=real(j)
            vrelaxn(:)=(s1*obc_vn(:,k)+s2*VVV(:,j,k))/real(spongen)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(spongen)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFKE_N(:,j,k)=amsue(:,0,k)*amsue(:,j,k)*lambda*(vrelaxn(:)-VVV(:,j,k))
          ENDDO
         ENDDO
       endif
      ENDIF

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
       	 DO k=1 ,ke
          DO j=je,je-sponges+1,-1
            s1=real(j-je+sponges)
            s2=real(je-j)
            urelaxs(:)=(s1*obc_us(:,k)+s2*UUU(:,j,k))/real(sponges)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(sponges)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            UFKO_S(:,j,k)=amsuo(:,je,k)*amsuo(:,j,k)*lambda*(urelaxs(:)-UUU(:,j,k))
          ENDDO

          DO j=je,je-sponges+1,-1
            s1=real(j-je+sponges)
            s2=real(je-j)
            vrelaxs(:)=(s1*obc_vs(:,k)+s2*VVV(:,j,k))/real(sponges)
            lambda=(s1*rel_uvb+s2*rel_uvi)/real(sponges)
            IF (lambda.ne.0.) THEN
              lambda=dt/lambda
            ELSE
              lambda=0.0
            ENDIF
            VFKE_S(:,j,k)=amsue(:,je,k)*amsue(:,j,k)*lambda*(vrelaxs(:)-VVV(:,j,k))
          ENDDO
         ENDDO
       endif
      endif

      UFKO=UFKO_W+UFKO_E+UFKO_N+UFKO_S
      VFKE=VFKE_W+VFKE_E+VFKE_N+VFKE_S

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then  
        DO i=0,spongew-1
          DO j=1,spongen
            if ((i+1)<j) then
              UFKO(i,j,:)=UFKO_W(i,j,:)
            else
              UFKO(i,j,:)=UFKO_N(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=0,spongen-1
            if (i<(j+1)) then
              VFKE(i,j,:)=VFKE_W(i,j,:)
            else
              VFKE(i,j,:)=VFKE_N(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then 
        DO i=0,spongew-1
          DO j=je,je-sponges+1,-1
            if ((i+1)<(je-j+1)) then
              UFKO(i,j,:)=UFKO_W(i,j,:)
            else
              UFKO(i,j,:)=UFKO_S(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              VFKE(i,j,:)=VFKE_W(i,j,:)
            else
              VFKE(i,j,:)=VFKE_S(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then 
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              UFKO(i,j,:)=UFKO_E(i,j,:)
            else
              UFKO(i,j,:)=UFKO_N(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=0,spongen-1
            if ((ie-i+1)<(j+1)) then
              VFKE(i,j,:)=VFKE_E(i,j,:)
            else
              VFKE(i,j,:)=VFKE_N(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then   
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              UFKO(i,j,:)=UFKO_E(i,j,:)
            else
              UFKO(i,j,:)=UFKO_S(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              VFKE(i,j,:)=VFKE_E(i,j,:)
            else
              VFKE(i,j,:)=VFKE_S(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      UUU0=UUU+UFKO
      VVV0=VVV+VFKE


      CALL bounds_exch(UUU0)
      CALL bounds_exch(VVV0)

     ! if ( East .AND. have_g_ie ) then
     !    write(IO_STDOUT,*) 'SPONGE_OBCS_UV WEST'
     !    write(IO_STDOUT,*) 'UUU0(ie,:,1)', UUU0(ie,:,1)
     !    write(IO_STDOUT,*) 'VVV0(ie,:,1)', VVV0(ie,:,1)
     !    write(IO_STDOUT,*) 'UUU0(ie-1,:,1)', UUU0(ie-1,:,1)
     !    write(IO_STDOUT,*) 'VVV0(ie-1,:,1)', VVV0(ie-1,:,1)         
     ! endif

     ! if ( West .AND. have_g_is ) then
     !    write(IO_STDOUT,*) 'SPONGE_OBCS_UV WEST'
     !    write(IO_STDOUT,*) 'UUU0(I_start,:,1)', UUU0(I_start,:,1)
     !    write(IO_STDOUT,*) 'VVV0(1,:,1)', VVV0(1,:,1)
     !    write(IO_STDOUT,*) 'UUU0(I_start+1,:,1)', UUU0(I_start+1,:,1)
     !    write(IO_STDOUT,*) 'VVV0(2,:,1)', VVV0(2,:,1)         
     ! endif


      END SUBROUTINE SPONGE_OBCS_UV


! H.H -- 2021.03.29
! relax OBC uv(3D) forcing to model fields by sponge layer
! with Exponent Function 
      SUBROUTINE SPONGE_OBCS_UV_EXP(UUU0,VVV0)

      REAL UUU0(I_start:IE,JE,KE),UUU(I_start:IE,JE,KE)
      REAL VVV0(IE,J_start:JE,KE),VVV(IE,J_start:JE,KE)
      
      REAL Lu_w(JE),Lu_e(JE),Lu_n(I_start:IE),Lu_s(I_start:IE)
      REAL Lv_w(J_start:JE),Lv_e(J_start:JE),Lv_n(IE),Lv_s(IE)  
         
      REAL Mu_w(JE),Mu_e(JE),Mu_n(I_start:IE),Mu_s(I_start:IE)
      REAL Mv_w(J_start:JE),Mv_e(J_start:JE),Mv_n(IE),Mv_s(IE)   
             
      UUU(:,:,:)=UUU0(:,:,:)
      VVV(:,:,:)=VVV0(:,:,:)
      
      UFKO=0.0
      UFKO_W=0.0
      UFKO_E=0.0
      UFKO_N=0.0
      UFKO_S=0.0

      VFKE=0.0
      VFKE_W=0.0
      VFKE_E=0.0
      VFKE_N=0.0
      VFKE_S=0.0
      
      Lu_w=0.0
      Lv_w=0.0
      Lu_e=0.0
      Lv_e=0.0
      Lu_n=0.0
      Lv_n=0.0
      Lu_s=0.0
      Lv_s=0.0
      
      Mu_w=0.0
      Mv_w=0.0
      Mu_e=0.0
      Mv_e=0.0
      Mu_n=0.0
      Mv_n=0.0
      Mu_s=0.0
      Mv_s=0.0      
      

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
      	  DO k=1,ke
           DO i=I_start,spongew-1

             Lu_w(:) = SUM(DLXU(I_start:i,:),DIM=1)            
             Mu_w(:) = exp(-4*Lu_w(:)/Lspu_w(:))/tao
             UFKO_W(i,:,k)=amsuo(I_start,:,k)*amsuo(i,:,k)*dt*Mu_w(:)*(obc_uw(:,k)-UUU(i,:,k))
           ENDDO
         

           DO i=1,spongew
             
             Lv_w(:) = SUM(DLXV(1:i,:),DIM=1)            
             Mv_w(:) = exp(-4*Lv_w(:)/Lspv_w(:))/tao
             
             VFKE_W(i,:,k)=amsue(1,:,k)*amsue(i,:,k)*dt*Mv_w(:)*(obc_vw(:,k)-VVV(i,:,k))
           ENDDO
          ENDDO
       endif
      endif

      

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
       	 DO k = 1,ke
          DO i=ie,ie-spongee+1,-1

            Lu_e(:) = SUM(DLXU(i:ie,:),DIM=1)            
            Mu_e(:) = exp(-4*Lu_e(:)/Lspu_e(:))/tao
                        
            UFKO_E(i,:,k)=amsuo(ie,:,k)*amsuo(i,:,k)*dt*Mu_e(:)*(obc_ue(:,k)-UUU(i,:,k))
          ENDDO

          DO i=ie,ie-spongee+1,-1

            Lv_e(:) = SUM(DLXV(i:ie,:),DIM=1)            
            Mv_e(:) = exp(-4*Lv_e(:)/Lspv_e(:))/tao

            VFKE_E(i,:,k)=amsue(ie,:,k)*amsue(i,:,k)*dt*Mv_e(:)*(obc_ve(:,k)-VVV(i,:,k))
          ENDDO
         endDO
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
        DO k = 1,ke 
          
          DO j=1,spongen
            Lu_n(:) = SUM(DLYU(:,1:j),DIM=2)            
            Mu_n(:) = exp(-4*Lu_n(:)/Lspu_n(:))/tao        
            UFKO_N(:,j,k)=amsuo(:,1,k)*amsuo(:,j,k)*dt*Mu_n(:)*(obc_un(:,k)-UUU(:,j,k))
          ENDDO

          DO j=J_start,spongen-1          
            Lv_n(:) = SUM(DLYV(:,J_start:j),DIM=2)            
            Mv_n(:) = exp(-4*Lv_n(:)/Lspv_n(:))/tao
            VFKE_N(:,j,k)=amsue(:,0,k)*amsue(:,j,k)*dt*Mv_n(:)*(obc_vn(:,k)-VVV(:,j,k))
          ENDDO
          
         ENDDO
       endif
      ENDIF

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
        DO k=1 ,ke
       
          DO j=je,je-sponges+1,-1       
            Lu_s(:) = SUM(DLYU(:,j:je),DIM=2)            
            Mu_s(:) = exp(-4*Lu_s(:)/Lspu_s(:))/tao
            UFKO_S(:,j,k)=amsuo(:,je,k)*amsuo(:,j,k)*dt*Mu_s(:)*(obc_us(:,k)-UUU(:,j,k))
          ENDDO

          DO j=je,je-sponges+1,-1
            Lv_s(:) = SUM(DLYV(:,j:je),DIM=2)            
            Mv_s(:) = exp(-4*Lv_s(:)/Lspv_s(:))/tao
            VFKE_S(:,j,k)=amsue(:,je,k)*amsue(:,j,k)*dt*Mv_s(:)*(obc_vs(:,k)-VVV(:,j,k))
          ENDDO
          
         ENDDO
       endif
      endif

      UFKO=UFKO_W+UFKO_E+UFKO_N+UFKO_S
      VFKE=VFKE_W+VFKE_E+VFKE_N+VFKE_S

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then  
        DO i=0,spongew-1
          DO j=1,spongen
            if ((i+1)<j) then
              UFKO(i,j,:)=UFKO_W(i,j,:)
            else
              UFKO(i,j,:)=UFKO_N(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=0,spongen-1
            if (i<(j+1)) then
              VFKE(i,j,:)=VFKE_W(i,j,:)
            else
              VFKE(i,j,:)=VFKE_N(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then 
        DO i=0,spongew-1
          DO j=je,je-sponges+1,-1
            if ((i+1)<(je-j+1)) then
              UFKO(i,j,:)=UFKO_W(i,j,:)
            else
              UFKO(i,j,:)=UFKO_S(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              VFKE(i,j,:)=VFKE_W(i,j,:)
            else
              VFKE(i,j,:)=VFKE_S(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then 
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              UFKO(i,j,:)=UFKO_E(i,j,:)
            else
              UFKO(i,j,:)=UFKO_N(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=0,spongen-1
            if ((ie-i+1)<(j+1)) then
              VFKE(i,j,:)=VFKE_E(i,j,:)
            else
              VFKE(i,j,:)=VFKE_N(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then   
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              UFKO(i,j,:)=UFKO_E(i,j,:)
            else
              UFKO(i,j,:)=UFKO_S(i,j,:)
            endif
          ENDDO
        ENDDO
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              VFKE(i,j,:)=VFKE_E(i,j,:)
            else
              VFKE(i,j,:)=VFKE_S(i,j,:)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      UUU0=UUU+UFKO
      VVV0=VVV+VFKE

      CALL bounds_exch(UUU0)
      CALL bounds_exch(VVV0)
      
      
    
      END SUBROUTINE SPONGE_OBCS_UV_EXP
      
! H.H -- 2021.03.29
! relax OBC T&S forcing to model fields by sponge layer
! with Exponent Function 
      SUBROUTINE SPONGE_OBCS_TS_EXP(TTT0,SSS0)
                 
      REAL TTT(IE,JE,KE),TTT0(IE,JE,KE),SSS(IE,JE,KE),SSS0(IE,JE,KE)
       
      REAL Lt_w(JE),Lt_e(JE),Lt_n(IE),Lt_s(IE)
      REAL Ls_w(JE),Ls_e(JE),Ls_n(IE),Ls_s(IE)  
         
      REAL Mt_w(JE),Mt_e(JE),Mt_n(IE),Mt_s(IE)
      REAL Ms_w(JE),Ms_e(JE),Ms_n(IE),Ms_s(IE)   
       
      TTT(:,:,:)=TTT0(:,:,:)
      SSS(:,:,:)=SSS0(:,:,:)      


      TFKO=0.0
      TFKO_W=0.0
      TFKO_E=0.0
      TFKO_N=0.0
      TFKO_S=0.0
      
      SFKO=0.0
      SFKO_W=0.0
      SFKO_E=0.0
      SFKO_N=0.0
      SFKO_S=0.0  
      
      Lt_w=0.0
      Ls_w=0.0
      Lt_e=0.0
      Ls_e=0.0
      Lt_n=0.0
      Ls_n=0.0
      Lt_s=0.0
      Ls_s=0.0   
      
      Mt_w=0.0
      Ms_w=0.0
      Mt_e=0.0
      Ms_e=0.0
      Mt_n=0.0
      Ms_n=0.0
      Mt_s=0.0
      Ms_s=0.0         

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
        DO k=1,ke      	
        DO i=1,spongew

          Lt_w(:) = SUM(DLXP(1:i,:),DIM=1)            
          Mt_w(:) = exp(-4*Lt_w(:)/Lspp_w(:))/tao
          TFKO_W(i,:,k)=weto(1,:,k)*weto(i,:,k)*dt*Mt_w(:)*(obc_tw(:,k)-TTT(i,:,k))
        
        ENDDO
        
        DO i=1,spongew
          
          Ls_w(:) = SUM(DLXP(1:i,:),DIM=1)            
          Ms_w(:) = exp(-4*Ls_w(:)/Lspp_w(:))/tao         
          SFKO_W(i,:,k)=weto(1,:,k)*weto(i,:,k)*dt*Ms_w(:)*(obc_sw(:,k)-SSS(i,:,k))
          
        ENDDO  
              
        ENDDO       
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
       
       DO k=1,ke
        DO i=ie,ie-spongee+1,-1
        
          Lt_e(:) = SUM(DLXP(i:ie,:),DIM=1)            
          Mt_e(:) = exp(-4*Lt_e(:)/Lspp_e(:))/tao
          TFKO_E(i,:,k)=weto(ie,:,k)*weto(i,:,k)*dt*Mt_e(:)*(obc_te(:,k)-TTT(i,:,k))
          
        ENDDO
 
        DO i=ie,ie-spongee+1,-1

          Ls_e(:) = SUM(DLXP(i:ie,:),DIM=1)            
          Ms_e(:) = exp(-4*Ls_e(:)/Lspp_e(:))/tao          
          SFKO_E(i,:,k)=weto(ie,:,k)*weto(i,:,k)*dt*Ms_e(:)*(obc_se(:,k)-SSS(i,:,k))
          
        ENDDO      
      
       ENDDO     
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
       DO k=1,ke
        DO j=1,spongen
        
           Lt_n(:) = SUM(DLYP(:,1:j),DIM=2)            
           Mt_n(:) = exp(-4*Lt_n(:)/Lspp_n(:))/tao  
           TFKO_N(:,j,k)=weto(:,1,k)*weto(:,j,k)*dt*Mt_n(:)*(obc_tn(:,k)-TTT(:,j,k))
           
        ENDDO

        DO j=1,spongen
        
           Ls_n(:) = SUM(DLYP(:,1:j),DIM=2)            
           Ms_n(:) = exp(-4*Ls_n(:)/Lspp_n(:))/tao  
           SFKO_N(:,j,k)=weto(:,1,k)*weto(:,j,k)*dt*Ms_n(:)*(obc_sn(:,k)-SSS(:,j,k))
           
        ENDDO
        
        ENDDO        
       endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
       DO k=1,ke 
             
        DO j=je,je-sponges+1,-1
        
           Lt_s(:) = SUM(DLYP(:,j:je),DIM=2)            
           Mt_s(:) = exp(-4*Lt_s(:)/Lspp_s(:))/tao
           TFKO_S(:,j,k)=weto(:,je,k)*weto(:,j,k)*dt*Mt_s(:)*(obc_ts(:,k)-TTT(:,j,k))
           
        ENDDO

        DO j=je,je-sponges+1,-1
        
           Ls_s(:) = SUM(DLYP(:,j:je),DIM=2)            
           Ms_s(:) = exp(-4*Lt_s(:)/Lspp_s(:))/tao
           SFKO_S(:,j,k)=weto(:,je,k)*weto(:,j,k)*dt*Ms_s(:)*(obc_ss(:,k)-SSS(:,j,k))
           
        ENDDO
        
        ENDDO        
        
       endif
      endif

      TFKO=TFKO_W+TFKO_E+TFKO_N+TFKO_S
      SFKO=SFKO_W+SFKO_E+SFKO_N+SFKO_S      
      

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then       	 
        DO i=1,spongew
          DO j=1,spongen
            if (i<j) then
              TFKO(i,j,:)=TFKO_W(i,j,:)
              SFKO(i,j,:)=SFKO_W(i,j,:)             
            else
              TFKO(i,j,:)=TFKO_N(i,j,:)
              SFKO(i,j,:)=SFKO_N(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then        	
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              TFKO(i,j,:)=TFKO_W(i,j,:)
              SFKO(i,j,:)=SFKO_W(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_S(i,j,:)
              SFKO(i,j,:)=SFKO_S(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then     
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              TFKO(i,j,:)=TFKO_E(i,j,:)
              SFKO(i,j,:)=SFKO_E(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_N(i,j,:)
              SFKO(i,j,:)=SFKO_N(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then        	
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              TFKO(i,j,:)=TFKO_E(i,j,:)
              SFKO(i,j,:)=SFKO_E(i,j,:)              
            else
              TFKO(i,j,:)=TFKO_S(i,j,:)
              SFKO(i,j,:)=SFKO_S(i,j,:)              
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      TTT0=TTT+TFKO
      SSS0=SSS+SFKO      

      CALL bounds_exch(TTT0)
      CALL bounds_exch(SSS0)

      END SUBROUTINE SPONGE_OBCS_TS_EXP    

! H.H -- 2021.03.30
! relax OBC Zeta forcing to model fields by sponge layer
! with Exponent Function 
! relax OBC eta forcing to model fields (with sponge layer)
      SUBROUTINE SPONGE_OBCS_ETA_EXP(ZZZ0)

      REAL ZZZ0(IE,JE),ZZZ(IE,JE)
      REAL Lp_w(JE),Lp_e(JE),Lp_n(IE),Lp_s(IE)         
      REAL Mp_w(JE),Mp_e(JE),Mp_n(IE),Mp_s(IE)
      
      ZZZ(:,:) = ZZZ0(:,:)

      Lp_w=0.0
      Lp_e=0.0
      Lp_n=0.0
      Lp_s=0.0  
      
      Mp_w=0.0
      Mp_e=0.0
      Mp_n=0.0
      Mp_s=0.0     

      ZFO=0.0
      ZFO_W=0.0
      ZFO_E=0.0
      ZFO_N=0.0
      ZFO_S=0.0

      ! west OB
      if ( West .AND. have_g_is) then       
       if (SpongeWest) then
        DO i=1,spongew
          Lp_w(:) = SUM(DLXP(1:i,:),DIM=1)            
          Mp_w(:) = exp(-4*Lp_w(:)/Lspp_w(:))/tao
          ZFO_W(i,:)=weto(1,:,1)*weto(i,:,1)*dt*Mp_w(:)*(obc_zw(:)-ZZZ(i,:))
        ENDDO
       endif
      endif

      ! east OB
      if ( East .AND. have_g_ie ) then      
       if (SpongeEast) then
        DO i=ie,ie-spongee+1,-1
          Lp_e(:) = SUM(DLXP(i:ie,:),DIM=1)            
          Mp_e(:) = exp(-4*Lp_e(:)/Lspp_e(:))/tao
          ZFO_E(i,:)=weto(ie,:,1)*weto(i,:,1)*dt*Mp_e(:)*(obc_ze(:)-ZZZ(i,:))
        ENDDO
       endif
      endif

      ! north OB
      if ( North .AND. have_g_js) then        
       if (SpongeNorth) then
        DO j=1,spongen
          Lp_n(:) = SUM(DLYP(:,1:j),DIM=2)            
          Mp_n(:) = exp(-4*Lp_n(:)/Lspp_n(:))/tao  
          ZFO_N(:,j)=weto(:,1,1)*weto(:,j,1)*dt*Mp_n(:)*(obc_zn(:)-ZZZ(:,j))
        ENDDO
       endif
      endif

      ! south OB
      if ( South .AND. have_g_je) then        
       if (SpongeSouth) then
        DO j=je,je-sponges+1,-1
          Lp_s(:) = SUM(DLYP(:,j:je),DIM=2)            
          Mp_s(:) = exp(-4*Lp_s(:)/Lspp_s(:))/tao
          ZFO_S(:,j)=weto(:,je,1)*weto(:,j,1)*dt*Mp_s(:)*(obc_zs(:)-ZZZ(:,j))
        ENDDO
       endif
      endif

      ZFO=ZFO_W+ZFO_E+ZFO_N+ZFO_S

      ! northwest OB
      if( North .and. West) then
       if (have_g_is .and. have_g_js) then 
       	if (SpongeWest .and. SpongeNorth) then       	 
        DO i=1,spongew
          DO j=1,spongen
            if (i<j) then
              ZFO(i,j)=ZFO_W(i,j)
            else
              ZFO(i,j)=ZFO_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southwest OB
      if( South .and. West) then      
       if (have_g_is .and. have_g_je) then
       	if (SpongeWest .and. SpongeSouth) then        	
        DO i=1,spongew
          DO j=je,je-sponges+1,-1
            if (i<(je-j+1)) then
              ZFO(i,j)=ZFO_W(i,j)
            else
              ZFO(i,j)=ZFO_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! northeast OB
      if( North .and. East) then         
       if (have_g_ie .and. have_g_js) then
       	if (SpongeEast .and. SpongeNorth) then     
        DO i=ie,ie-spongee+1,-1
          DO j=1,spongen
            if ((ie-i+1)<j) then
              ZFO(i,j)=ZFO_E(i,j)
            else
              ZFO(i,j)=ZFO_N(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ! southeast OB
      if( South .and. East) then        
       if (have_g_ie .and. have_g_je) then
        if (SpongeEast .and. SpongeSouth) then        	
        DO i=ie,ie-spongee+1,-1
          DO j=je,je-sponges+1,-1
            if ((ie-i+1)<(je-j+1)) then
              ZFO(i,j)=ZFO_E(i,j)
            else
              ZFO(i,j)=ZFO_S(i,j)
            endif
          ENDDO
        ENDDO
        endif
       endif
      endif

      ZZZ0=ZZZ+ZFO

      CALL bounds_exch(ZZZ0)

      
       
      END SUBROUTINE SPONGE_OBCS_ETA_EXP

! check if OBC data has been changed in processes other than OBC module
      SUBROUTINE CHECK_OBCS_ETA(i_step)

      ! west OB
      if ( West ) then      
      if (have_g_is) then
        do j=1,je
          if (abs(obc_zw(j)-zo(1,j))>obcerr)                            &
     &                              write(IO_STDOUT,1),'ZW','j',j,i_step
        enddo
      endif
      endif

      ! east OB
      if ( East ) then         
      if (have_g_ie) then
        do j=1,je
          if (abs(obc_ze(j)-zo(ie,j))>obcerr)                           &
     &                              write(IO_STDOUT,1),'ZE','j',j,i_step
        enddo
      endif
      endif
      
      ! north OB
      if ( North ) then         
      if (have_g_js) then
        do i=1,ie
          if (abs(obc_zn(i)-zo(i,1))>obcerr)                            &
     &                              write(IO_STDOUT,1),'ZN','i',i,i_step
        enddo
      endif
      endif

      ! south OB
      if ( South ) then         
      if (have_g_je) then
        do i=1,ie
          if (abs(obc_zs(i)-zo(i,je))>obcerr)                           &
     &                              write(IO_STDOUT,1),'ZS','i',i,i_step
        enddo
      endif
      endif

1     FORMAT('CHECK ETA OBC_',A2,': ',A1,'=',I3,',t=',I8)

      END SUBROUTINE CHECK_OBCS_ETA


! check if OBC data has been changed in processes other than OBC module
      SUBROUTINE CHECK_OBCS_UV(i_step)

      ! west OB
      if ( West ) then       
       if (have_g_is) then
        do k=1,ke
          do j=J_start,je
            if (abs(obc_vw(j,k)-voe(1,j,k))>obcerr)                     &
     &                            write(IO_STDOUT,3),'VW','j',j,k,i_step
          enddo
          do j=1,je
            if (abs(obc_uw(j,k)-uoo(0,j,k))>obcerr)                     &
     &                            write(IO_STDOUT,3),'UW','j',j,k,i_step
          enddo
        enddo
       endif
      endif

      ! east OB
      if ( East ) then       
       if (have_g_ie) then
        do k=1,ke
          do j=J_start,je
            if (abs(obc_ve(j,k)-voe(ie,j,k))>obcerr)                    &
     &                            write(IO_STDOUT,3),'VE','j',j,k,i_step
          enddo
          do j=1,je
            if (abs(obc_ue(j,k)-uoo(ie,j,k))>obcerr)                    &
     &                            write(IO_STDOUT,3),'UE','j',j,k,i_step
          enddo
        enddo
       endif
      endif

      ! north OB
      if ( North ) then       
       if (have_g_js) then
        do k=1,ke
          do i=I_start,ie
            if (abs(obc_un(i,k)-uoo(i,1,k))>obcerr)                     &
     &                            write(IO_STDOUT,3),'UN','i',i,k,i_step
          enddo
          do i=1,ie
            if (abs(obc_vn(i,k)-voe(i,0,k))>obcerr)                     &
     &                            write(IO_STDOUT,3),'VN','i',i,k,i_step
          enddo
        enddo
       endif
      endif

      ! south OB
      if ( South ) then       
       if (have_g_je) then
        do k=1,ke
          do i=I_start,ie
            if (abs(obc_us(i,k)-uoo(i,je,k))>obcerr)                    &
     &                            write(IO_STDOUT,3),'US','i',i,k,i_step
          enddo
          do i=1,ie
            if (abs(obc_vs(i,k)-voe(i,je,k))>obcerr)                    &
     &                            write(IO_STDOUT,3),'VS','i',i,k,i_step
          enddo
        enddo
       endif
      endif

3     FORMAT('CHECK UV OBC_',A2,': ',A1,'=',I3,',k=',I2,',t=',I8)

      END SUBROUTINE CHECK_OBCS_UV


! check if OBC data has been changed in processes other than OBC module
      SUBROUTINE CHECK_OBCS_TS(i_step)

      ! west OB
      if ( West ) then       
       if (have_g_is) then
        do k=1,ke
          do j=1,je
            if (abs(obc_tw(j,k)-tho(1,j,k))>obcerr)                     &
     &                            write(IO_STDOUT,2),'TW','j',j,k,i_step
            if (abs(obc_sw(j,k)-sao(1,j,k))>obcerr)                     &
     &                            write(IO_STDOUT,2),'SW','j',j,k,i_step
          enddo
        enddo
       endif
      endif

      ! east OB
      if ( East ) then       
       if (have_g_ie) then
        do k=1,ke
          do j=1,je
            if (abs(obc_te(j,k)-tho(ie,j,k))>obcerr)                    &
     &                            write(IO_STDOUT,2),'TE','j',j,k,i_step
            if (abs(obc_se(j,k)-sao(ie,j,k))>obcerr)                    &
     &                            write(IO_STDOUT,2),'SE','j',j,k,i_step
          enddo
        enddo
       endif
      endif

      ! north OB
      if ( North ) then       
       if (have_g_js) then
        do k=1,ke
          do i=1,ie
            if (abs(obc_tn(i,k)-tho(i,1,k))>obcerr)                     &
     &                            write(IO_STDOUT,2),'TN','i',i,k,i_step
            if (abs(obc_sn(i,k)-sao(i,1,k))>obcerr)                     &
     &                            write(IO_STDOUT,2),'SN','i',i,k,i_step
          enddo
        enddo
       endif
      endif

      ! south OB
      if ( South ) then       
       if (have_g_je) then
        do k=1,ke
          do i=1,ie
            if (abs(obc_ts(i,k)-tho(i,je,k))>obcerr)                    &
     &                            write(IO_STDOUT,2),'TS','i',i,k,i_step
            if (abs(obc_ss(i,k)-sao(i,je,k))>obcerr)                    &
     &                            write(IO_STDOUT,2),'SS','i',i,k,i_step
          enddo
        enddo
       endif
      endif

2     FORMAT('CHECK TS OBC_',A2,': ',A1,'=',I3,',k=',I2,',t=',I8)

      END SUBROUTINE CHECK_OBCS_TS


      END MODULE MO_OBCS


